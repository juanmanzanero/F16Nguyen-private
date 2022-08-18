#include "lion/thirdparty/include/cppad/ipopt/solve.hpp"
#include "lion/math/ipopt_cppad_handler.hpp"

template<typename Dynamic_model_t>
inline Optimal_laptime<Dynamic_model_t>::Optimal_laptime(const std::vector<scalar>& s_, const bool is_closed_, const bool is_direct_,
    const Dynamic_model_t& car, const std::vector<std::array<scalar,Dynamic_model_t::NSTATE>>& q0, 
    const std::vector<std::array<scalar,Dynamic_model_t::NALGEBRAIC>>& qa0,
    Control_variables<> control_variables_0, 
    const Options opts)
: options(opts), is_closed(is_closed_), is_direct(is_direct_), warm_start(false), 
  s(s_), q(q0), qa(qa0), control_variables(control_variables_0.check()),
  optimization_data()
{
    // (1) Check inputs
    check_inputs(car);

    // (4) Compute
    compute(car);
}


template<typename Dynamic_model_t>
inline Optimal_laptime<Dynamic_model_t>::Optimal_laptime(const std::vector<scalar>& s_, const bool is_closed_, const bool is_direct_,
    const Dynamic_model_t& car, const std::vector<std::array<scalar,Dynamic_model_t::NSTATE>>& q0, 
    const std::vector<std::array<scalar,Dynamic_model_t::NALGEBRAIC>>& qa0,
    Control_variables<> control_variables_0,
    const std::vector<scalar>& zl,
    const std::vector<scalar>& zu,
    const std::vector<scalar>& lambda,
    const Options opts)
: options(opts), is_closed(is_closed_), is_direct(is_direct_), warm_start(true), 
  s(s_), q(q0), qa(qa0), control_variables(control_variables_0.check()),
  optimization_data({.zl{zl},.zu{zu},.lambda{lambda}})
{
    // (1) Check inputs
    check_inputs(car);
      
    // (5) Compute
    compute(car);
}


template<typename Dynamic_model_t>
void Optimal_laptime<Dynamic_model_t>::check_inputs(const Dynamic_model_t& car)
{
    if ( s.size() <= 1 )
        throw std::runtime_error("Provide at least two values of arclength");

    // (1) Compute number of elements and points
    n_points = s.size();
    n_elements = (is_closed ? n_points : n_points - 1);

    // (2) Verify the vector or arclength
    const scalar& L = car.get_road().track_length();

    if (is_closed)
    {
        if (std::abs(s.front()) > 1.0e-12)
            throw std::runtime_error("In closed circuits, s[0] should be 0.0");

        if (s.back() > L - 1.0e-10)
            throw std::runtime_error("In closed circuits, s[end] should be < track_length");
    }
    else
    {
        if (s[0] < -1.0e-12)
            throw std::runtime_error("s[0] must be >= 0");

        if (s.back() > L)
            throw std::runtime_error("s[end] must be <= L");
    }

    // If closed, replace the initial arclength by 0
    if (is_closed)
    {   
        s.front() = 0.0;
    }

    // (3) Set the initial condition

    // (3.1) State
    if ( q.size() != n_points )
        throw std::runtime_error("q must have size of n_points");

    // (3.2) Algebraic state
    if ( qa.size() != n_points )
        throw std::runtime_error("qa must have size of n_points");
    
    // (3.3) Controls: check size only for full-mesh controls
    for (const auto& control_variable : control_variables )
    {
        if ( control_variable.optimal_control_type == FULL_MESH )
        {
            if ( control_variable.u.size() != n_points )
                throw std::runtime_error("control_variables.u must have the size of n_points");            

            if ( is_direct )
            {
                if ( control_variable.dudt.size() != 0 )
                    throw std::runtime_error("In direct simulations, control_variables.dudt must be empty");            
            }
            else
            {   
                if ( control_variable.dudt.size() != n_points )
                    throw std::runtime_error("control_variables.dudt must have the size of n_points");            
            }
        }
    }
}



template<typename Dynamic_model_t>
inline Optimal_laptime<Dynamic_model_t>::Optimal_laptime(Xml_document& doc)
{
    Xml_element root = doc.get_root_element();

    if ( root.get_attribute("type") == "closed" )
        is_closed = true;
    else if ( root.get_attribute("type") == "open" )
        is_closed = false;
    else
        throw std::runtime_error("Incorrect track type, should be \"open\" or \"closed\"");

    if ( root.get_attribute("is_direct") == "true" )
        is_direct = true;
    else if ( root.get_attribute("is_direct") == "false" )
        is_direct = false;
    else
        throw std::runtime_error("Incorrect \"is_direct\" attribute, should be \"true\" or \"false\"");

    const auto [key_name, q_names, qa_names, u_names] = Dynamic_model_t::get_state_and_control_names();

    // Get the data
    n_points = std::stoi(root.get_attribute("n_points"));
    n_elements = (is_closed ? n_points : n_points - 1);

    laptime = root.get_child("laptime").get_value(scalar());

    s = root.get_child("arclength").get_value(std::vector<scalar>());

    // Get state
    q = std::vector<std::array<scalar,Dynamic_model_t::NSTATE>>(n_points);
    for (size_t i = 0; i < Dynamic_model_t::NSTATE; ++i)
    {
        std::vector<scalar> data_in = root.get_child(q_names[i]).get_value(std::vector<scalar>());
        for (size_t j = 0; j < n_points; ++j)
            q[j][i] = data_in[j];
    }

    // Get algebraic states
    qa = std::vector<std::array<scalar,Dynamic_model_t::NALGEBRAIC>>(n_points);
    for (size_t i = 0; i < Dynamic_model_t::NALGEBRAIC; ++i)
    {
        std::vector<scalar> data_in = root.get_child(qa_names[i]).get_value(std::vector<scalar>());
        for (size_t j = 0; j < n_points; ++j)
            qa[j][i] = data_in[j];
    }

    // Get controls
    control_variables = Control_variables<>{};
    for (size_t i = 0; i < Dynamic_model_t::NCONTROL; ++i)
    {
        auto control_var_element = root.get_child("control_variables/" + u_names[i]);
        
        // Get optimal control type attribute
        auto optimal_control_type_str = control_var_element.get_attribute("optimal_control_type");

        if ( optimal_control_type_str == "dont optimize" )
            control_variables[i].optimal_control_type = DONT_OPTIMIZE;
        else if ( optimal_control_type_str == "constant" )
            control_variables[i].optimal_control_type = CONSTANT;
        else if ( optimal_control_type_str == "hypermesh" )
            control_variables[i].optimal_control_type = HYPERMESH;
        else if ( optimal_control_type_str == "full-mesh" )
            control_variables[i].optimal_control_type = FULL_MESH;
        else
            throw std::runtime_error("optimal_control_type attribute not recognize. Options are: 'dont optimize', 'constant', 'hypermesh', 'full-mesh'");

        // Get value
        control_variables[i].u = root.get_child("control_variables/" + u_names[i] + "/values").get_value(std::vector<scalar>());

        // Get hypermesh
        if ( control_variables[i].optimal_control_type == HYPERMESH )
            control_variables[i].s_hypermesh = root.get_child("control_variables/" + u_names[i] + "/hypermesh").get_value(std::vector<scalar>());

        // Get derivative value if present
        if ( !is_direct && root.has_child("control_variables/" + u_names[i] + "/derivatives") )
            control_variables[i].dudt = root.get_child("control_variables/" + u_names[i] + "/derivatives").get_value(std::vector<scalar>());

        // Get hypermesh 
        if ( control_variables[i].optimal_control_type == HYPERMESH )
            control_variables[i].s_hypermesh = root.get_child("control_variables/" + u_names[i] + "/hypermesh").get_value(std::vector<scalar>());
    }

    // check them and compute its statistics
    control_variables.check();

    // Get x
    x_coord = root.get_child("x").get_value(std::vector<scalar>());
    y_coord = root.get_child("y").get_value(std::vector<scalar>());
    psi     = root.get_child("psi").get_value(std::vector<scalar>());

    // Get optimization data
    optimization_data.zl     = root.get_child("optimization_data").get_child("zl").get_value(std::vector<scalar>());
    optimization_data.zu     = root.get_child("optimization_data").get_child("zu").get_value(std::vector<scalar>());
    optimization_data.lambda = root.get_child("optimization_data").get_child("lambda").get_value(std::vector<scalar>());
}


template<typename Dynamic_model_t>
inline void Optimal_laptime<Dynamic_model_t>::compute(const Dynamic_model_t& car)
{
    if ( is_direct )
    {
        if ( is_closed )
            compute_direct<true>(car);
        else
            compute_direct<false>(car);
    }
    else
    {
        if ( is_closed )
            compute_derivative<true>(car);
        else
            compute_derivative<false>(car);
    }
}


template<typename Dynamic_model_t>
template<bool isClosed>
inline void Optimal_laptime<Dynamic_model_t>::compute_direct(const Dynamic_model_t& car) 
{
    class Counter
    {
     public:
        Counter() : k(0) {}

        Counter& increment() { ++k; return *this;}
    
        operator size_t() const { return k; }

     private:
        size_t k;
    };

    // (1) Get variable bounds and default control variables
    const auto [__ignore, q_lb, q_ub, ___ignore, qa_lb, qa_ub, ____ignore, u_lb, u_ub] = car.get_state_and_control_upper_lower_and_default_values();
    (void) __ignore;
    (void) ___ignore;
    (void) ____ignore;

    // (1) Construct starting vector of control variables
    auto u0 = control_variables.control_array_at_s(car, 0, s.front());

    // (2) Construct fitness function functor
    FG_direct<isClosed> fg(n_elements,n_points,car,s,q.front(),qa.front(),u0,control_variables, options.sigma);

    // (3) Construct vectors of initial optimization point, and variable upper/lower bounds
    std::vector<scalar> x0(fg.get_n_variables(),0.0);
    std::vector<scalar> x_lb(fg.get_n_variables(), std::numeric_limits<scalar>::lowest());
    std::vector<scalar> x_ub(fg.get_n_variables(), std::numeric_limits<scalar>::max());

    Counter k;
    constexpr const size_t offset = isClosed ? 0 : 1;

    for (size_t i = offset; i < n_points; ++i)
    {
        // (3.1) Set state and initial condition from start to ITIME
        for (size_t j = 0; j < Dynamic_model_t::Road_type::ITIME; ++j)    
        {
            x0[k] = q[i][j];
            x_lb[k] = q_lb[j];
            x_ub[k] = q_ub[j];
            k.increment();
        }

        // (3.2) Set state to IN. Assert that ITIME = IN - 1
        assert(Dynamic_model_t::Road_type::ITIME == ( Dynamic_model_t::Road_type::IN - 1 ) );

        x0[k] = q[i][Dynamic_model_t::Road_type::IN];
        x_lb[k] = -car.get_road().get_left_track_limit(s[i]);
        x_ub[k] =  car.get_road().get_right_track_limit(s[i]);
        k.increment();

        // (3.3) Set state after IN
        for (size_t j = Dynamic_model_t::Road_type::IN+1; j < Dynamic_model_t::NSTATE; ++j)    
        {
            x0[k] = q[i][j];
            x_lb[k] = q_lb[j];
            x_ub[k] = q_ub[j];
            k.increment();
        }

        // (3.4) Set algebraic variables
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
        {
            x0[k] = qa[i][j];
            x_lb[k] = qa_lb[j];
            x_ub[k] = qa_ub[j];
            k.increment();
        }

        // (3.5) Set full mesh control variables
        for (size_t j = 0; j < Dynamic_model_t::NCONTROL; ++j)
        {
            if ( control_variables[j].optimal_control_type == FULL_MESH )
            {
                x0[k] = control_variables[j].u[i];
                x_lb[k] = u_lb[j];
                x_ub[k] = u_ub[j];
                k.increment();
            }
        }
    }

    // (3.6) Add the control variables corresponding to constant optimizations and hypermesh optimizations
    for (size_t j = 0; j < Dynamic_model_t::NCONTROL; ++j)
    {
        switch (control_variables[j].optimal_control_type)
        {
         case (CONSTANT):
            for (size_t i = 0; i < control_variables[j].u.size(); ++i)
            {
                x0[k] = control_variables[j].u[i];
                x_lb[k] = u_lb[j];
                x_ub[k] = u_ub[j];
                k.increment();
            } 
            break;

         case (HYPERMESH):
            for (size_t i = 0; i < control_variables[j].u.size(); ++i)
            {
                x0[k] = control_variables[j].u[i];
                x_lb[k] = u_lb[j];
                x_ub[k] = u_ub[j];
                k.increment();
            }    
            break;
         default:
            break;
        }
    }

    // (4) Construct vector of upper/lower bounds for constraints
    std::vector<scalar> c_lb(fg.get_n_constraints(),0.0);
    std::vector<scalar> c_ub(fg.get_n_constraints(),0.0);
    Counter kc;
    for (size_t i = 1; i < n_points; ++i)
    {
        // Equality constraints: --------------- 

        // (4.1) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (before time)
        for (size_t j = 0; j < Dynamic_model_t::Road_type::ITIME; ++j)
        {
            c_lb[kc] = 0.0;
            c_ub[kc] = 0.0;
            kc.increment();
        }

        // (4.2) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (after time)
        for (size_t j = Dynamic_model_t::Road_type::ITIME+1; j < Dynamic_model_t::NSTATE; ++j)
        {
            c_lb[kc] = 0.0;
            c_ub[kc] = 0.0;
            kc.increment();
        }

        // (4.3) Algebraic constraints: dqa^{i} = 0.0
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
        {
            c_lb[kc] = 0.0;
            c_ub[kc] = 0.0;
            kc.increment();
        }
        
        // (4.4) Extra constraints per point
        auto [c_extra_lb, c_extra_ub] = car.optimal_laptime_extra_constraints_bounds(s[i]); 
        for (size_t j = 0; j < Dynamic_model_t::N_OL_EXTRA_CONSTRAINTS; ++j)
        {
            c_lb[kc] = c_extra_lb[j];
            c_ub[kc] = c_extra_ub[j];
            kc.increment();
        }
    }

    // (5) Add the periodic element if track is closed
    if constexpr (isClosed)
    {
        // Equality constraints: --------------- 

        // (5.1) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (before time)
        for (size_t j = 0; j < Dynamic_model_t::Road_type::ITIME; ++j)
        {
            c_lb[kc] = 0.0;
            c_ub[kc] = 0.0;
            kc.increment();
        }

        // (5.2) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (after time)
        for (size_t j = Dynamic_model_t::Road_type::ITIME+1; j < Dynamic_model_t::NSTATE; ++j)
        {
            c_lb[kc] = 0.0;
            c_ub[kc] = 0.0;
            kc.increment();
        }

        // (5.3) Algebraic constraints: dqa^{i} = 0.0
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
        {
            c_lb[kc] = 0.0;
            c_ub[kc] = 0.0;
            kc.increment();
        }
        
        // (5.4) Extra constraints per point
        auto [c_extra_lb, c_extra_ub] = car.optimal_laptime_extra_constraints_bounds(s[0]); 
        for (size_t j = 0; j < Dynamic_model_t::N_OL_EXTRA_CONSTRAINTS; ++j)
        {
            c_lb[kc] = c_extra_lb[j];
            c_ub[kc] = c_extra_ub[j];
            kc.increment();
        }
    }

    // (6) Check the dimensions of the counter used to fill the vectors
    if ( k != fg.get_n_variables() )
    {
        std::ostringstream s_out;
        s_out << "k(=" << k << ") != fg.get_n_variables()(=" << fg.get_n_variables() << ")" ;
        throw std::runtime_error(s_out.str());
    }
    
    if ( kc != fg.get_n_constraints() )
    {
        std::ostringstream s_out;
        s_out << "kc(=" << kc << ") != fg.get_n_constraints()(=" << fg.get_n_constraints() << ")" ;
        throw std::runtime_error(s_out.str());
    }

    // (7) Warm-start: simply check the dimensions
    if ( warm_start )
    {
        if ( optimization_data.zl.size() != fg.get_n_variables() )
        {
            std::ostringstream s_out;
            s_out << "size of zl should be " << fg.get_n_variables() << " but it is " << optimization_data.zl.size();
            throw std::runtime_error(s_out.str());
        }

        if ( optimization_data.zu.size() != fg.get_n_variables() )
        {
            std::ostringstream s_out;
            s_out << "size of zu should be " << fg.get_n_variables() << " but it is " << optimization_data.zu.size();
            throw std::runtime_error(s_out.str());
        }

        if ( optimization_data.lambda.size() != fg.get_n_constraints() )
        {
            std::ostringstream s_out;
            s_out << "size of lambda should be " << fg.get_n_constraints() << " but it is " << optimization_data.lambda.size();
            throw std::runtime_error(s_out.str());
        }
    } 

    // (8) Run optimization

    // (8.1) Prepare options
    std::string ipoptoptions;
    ipoptoptions += "Integer print_level  ";
    ipoptoptions += std::to_string(options.print_level);
    ipoptoptions += "\n";
    ipoptoptions += "Integer max_iter ";
    ipoptoptions += std::to_string(options.maximum_iterations);
    ipoptoptions += "\n";
    ipoptoptions += "String  sb           yes\n";
    ipoptoptions += "String nlp_scaling_method gradient-based\n";
    ipoptoptions += "Numeric nlp_scaling_max_gradient 1.0\n";
    ipoptoptions += "Sparse true forward\n";
    ipoptoptions += "Retape true\n";
    ipoptoptions += "Numeric tol          1e-06\n";
    ipoptoptions += "Numeric constr_viol_tol  1e-06\n";
    ipoptoptions += "Numeric acceptable_tol  1e-4\n";

    // (8.2) Return object
    CppAD::ipopt_cppad_result<std::vector<scalar>> result;

    // (8.3) Solve the problem
    if ( !warm_start )
        CppAD::ipopt_cppad_solve<std::vector<scalar>, FG_direct<isClosed>>(ipoptoptions, x0, x_lb, x_ub, c_lb, c_ub, fg, result);
    else
        CppAD::ipopt_cppad_solve<std::vector<scalar>, FG_direct<isClosed>>(ipoptoptions, x0, x_lb, x_ub, c_lb, c_ub, 
            optimization_data.lambda, optimization_data.zl, optimization_data.zu, fg, result);
    

    // (8.4) Check success flag
    success = result.status == CppAD::ipopt_cppad_result<std::vector<scalar>>::success; 
    iter_count = result.iter_count;

    if ( !success && options.throw_if_fail )
    {
        throw std::runtime_error("Optimization did not succeed");
    }

    // (9) Export the solution

    // (9.1) Load the latest solution vector to the fitness function object
    std::vector<CppAD::AD<scalar>> x_final(result.x.size(),0.0),fg_final(fg.get_n_constraints()+1,0.0);
    std::copy(result.x.cbegin(), result.x.cend(), x_final.begin());
    fg(fg_final,x_final);

    assert(fg.get_states().size()           == n_points);
    assert(fg.get_algebraic_states().size() == n_points);

    // (9.2) Export states
    q = std::vector<std::array<scalar,Dynamic_model_t::NSTATE>>(n_points);
    
    for (size_t i = 0; i < n_points; ++i)
        for (size_t j = 0; j < Dynamic_model_t::NSTATE; ++j)
            q[i][j] = Value(fg.get_state(i)[j]);

    // (9.3) Export algebraic states
    qa = std::vector<std::array<scalar,Dynamic_model_t::NALGEBRAIC>>(n_points);
    
    for (size_t i = 0; i < n_points; ++i)
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
            qa[i][j] = Value(fg.get_algebraic_state(i)[j]);

    // (9.4) Export control variables
    control_variables.clear();
    
    for (size_t j = 0; j < Dynamic_model_t::NCONTROL; ++j)
    {
        switch (control_variables[j].optimal_control_type) 
        {
         case (CONSTANT):
            std::transform(fg.get_control(j).u.cbegin(), fg.get_control(j).u.cend(), control_variables[j].u.begin(), 
                [](const auto& u_cppad) -> auto { return Value(u_cppad); });
            break;
        
         case (HYPERMESH):
            std::transform(fg.get_control(j).u.cbegin(), fg.get_control(j).u.cend(), control_variables[j].u.begin(), 
                [](const auto& u_cppad) -> auto { return Value(u_cppad); });
            break;

         case (FULL_MESH):
            std::transform(fg.get_control(j).u.cbegin(), fg.get_control(j).u.cend(), control_variables[j].u.begin(), 
                [](const auto& u_cppad) -> auto { return Value(u_cppad); });
            break;

         default:
            // do nothing
            break;
        }
    }

    // (9.5) Export time, x, y, and psi
    const scalar& L = car.get_road().track_length();

    x_coord = std::vector<scalar>(fg.get_states().size());
    y_coord = std::vector<scalar>(fg.get_states().size());
    psi = std::vector<scalar>(fg.get_states().size());

    // (9.5.1) Compute the first point
    auto u_i = fg.get_controls().control_array_at_s(car, 0, s.front());
    auto dtimeds_first = fg.get_car()(fg.get_state(0),fg.get_algebraic_state(0),u_i,s[0]).first[Dynamic_model_t::Road_type::ITIME];
    auto dtimeds_prev = dtimeds_first;

    x_coord.front() = Value(fg.get_car().get_road().get_x());
    y_coord.front() = Value(fg.get_car().get_road().get_y());
    psi.front() = Value(fg.get_car().get_road().get_psi());

    // (9.5.2) Compute the rest of the points
    for (size_t i = 1; i < fg.get_states().size(); ++i)
    {
        u_i = fg.get_controls().control_array_at_s(car, i, s[i]);
        const auto dtimeds = fg.get_car()(fg.get_state(i),fg.get_algebraic_state(i),u_i,s[i]).first[Dynamic_model_t::Road_type::ITIME];
        q[i][Dynamic_model_t::Road_type::ITIME] = q[i-1][Dynamic_model_t::Road_type::ITIME] 
            + Value((s[i]-s[i-1])*(options.sigma*dtimeds + (1.0-options.sigma)*dtimeds_prev));
        dtimeds_prev = dtimeds;

        x_coord[i] = Value(fg.get_car().get_road().get_x());
        y_coord[i] = Value(fg.get_car().get_road().get_y());
        psi[i]     = Value(fg.get_car().get_road().get_psi());
    }

    // (9.5.4) Compute the laptime
    laptime = q.back()[Dynamic_model_t::Road_type::ITIME];

    if (isClosed)
        laptime += Value((L-s.back())*((1.0-options.sigma)*dtimeds_prev + options.sigma*dtimeds_first));

    // (9.6) Save optimization data
    optimization_data.zl     = result.zl;
    optimization_data.zu     = result.zu;
    optimization_data.lambda = result.lambda;
}

template<typename Dynamic_model_t>
template<bool isClosed>
inline void Optimal_laptime<Dynamic_model_t>::compute_derivative(const Dynamic_model_t& car) 
{
    class Counter
    {
     public:
        Counter() : k(0) {}

        Counter& increment() { ++k; return *this;}
    
        operator size_t() const { return k; }

     private:
        size_t k;
    };

    // (1) Get variable bounds and default control variables
    const auto [__ignore, q_lb, q_ub, ___ignore, qa_lb, qa_ub, ____ignore, u_lb, u_ub] = car.get_state_and_control_upper_lower_and_default_values();
    (void) __ignore;
    (void) ___ignore;
    (void) ____ignore;

    const auto [dudt_lb, dudt_ub] = car.optimal_laptime_derivative_control_bounds();

    // (1) Construct starting vector of control variables
    auto [u0, dudt0] = control_variables.control_array_and_derivative_at_s(car, 0, s.front());

    // (2) Construct fitness function functor
    FG_derivative<isClosed> fg(n_elements,n_points,car,s,q.front(),qa.front(),u0,dudt0,control_variables,options.sigma);

    // (3) Construct vectors of initial optimization point, and variable upper/lower bounds
    std::vector<scalar> x0(fg.get_n_variables(),0.0);
    std::vector<scalar> x_lb(fg.get_n_variables(), std::numeric_limits<scalar>::lowest());
    std::vector<scalar> x_ub(fg.get_n_variables(), std::numeric_limits<scalar>::max());

    Counter k;
    constexpr const size_t offset = isClosed ? 0 : 1;

    for (size_t i = offset; i < n_points; ++i)
    {
        // (3.1) Set state and initial condition from start to ITIME
        for (size_t j = 0; j < Dynamic_model_t::Road_type::ITIME; ++j)    
        {
            x0[k] = q[i][j];
            x_lb[k] = q_lb[j];
            x_ub[k] = q_ub[j];
            k.increment();
        }

        // (3.2) Set state to IN. Assert that ITIME = IN - 1
        assert(Dynamic_model_t::Road_type::ITIME == ( Dynamic_model_t::Road_type::IN - 1 ) );

        x0[k] = q[i][Dynamic_model_t::Road_type::IN];
        x_lb[k] = -car.get_road().get_left_track_limit(s[i]);
        x_ub[k] =  car.get_road().get_right_track_limit(s[i]);
        k.increment();

        // (3.3) Set state after IN
        for (size_t j = Dynamic_model_t::Road_type::IN+1; j < Dynamic_model_t::NSTATE; ++j)    
        {
            x0[k] = q[i][j];
            x_lb[k] = q_lb[j];
            x_ub[k] = q_ub[j];
            k.increment();
        }

        // (3.4) Set algebraic variables
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
        {
            x0[k] = qa[i][j];
            x_lb[k] = qa_lb[j];
            x_ub[k] = qa_ub[j];
            k.increment();
        }

        // (3.5) Set full mesh control variables
        for (size_t j = 0; j < Dynamic_model_t::NCONTROL; ++j)
        {
            if ( control_variables[j].optimal_control_type == FULL_MESH )
            {
                x0[k] = control_variables[j].u[i];
                x_lb[k] = u_lb[j];
                x_ub[k] = u_ub[j];
                k.increment();
            }
        }

        // (3.6) Set full mesh control variables derivatives
        for (size_t j = 0; j < Dynamic_model_t::NCONTROL; ++j)
        {
            if ( control_variables[j].optimal_control_type == FULL_MESH )
            {
                x0[k] = control_variables[j].dudt[i];
                x_lb[k] = dudt_lb[j];
                x_ub[k] = dudt_ub[j];
                k.increment();
            }
        }
    }

    // (3.7) Add the control variables corresponding to constant optimizations and hypermesh optimizations
    for (size_t j = 0; j < Dynamic_model_t::NCONTROL; ++j)
    {
        switch (control_variables[j].optimal_control_type)
        {
         case (CONSTANT):
            for (size_t i = 0; i < control_variables[j].u.size(); ++i)
            {
                x0[k] = control_variables[j].u[i];
                x_lb[k] = u_lb[j];
                x_ub[k] = u_ub[j];
                k.increment();
            } 

         case (HYPERMESH):
            for (size_t i = 0; i < control_variables[j].u.size(); ++i)
            {
                x0[k] = control_variables[j].u[i];
                x_lb[k] = u_lb[j];
                x_ub[k] = u_ub[j];
                k.increment();
            }    
         default:
            break;
        }
    }

    // (4) Construct vector of upper/lower bounds for constraints
    std::vector<scalar> c_lb(fg.get_n_constraints(),0.0);
    std::vector<scalar> c_ub(fg.get_n_constraints(),0.0);
    Counter kc;
    for (size_t i = 1; i < n_points; ++i)
    {
        // Equality constraints: --------------- 

        // (4.1) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (before time)
        for (size_t j = 0; j < Dynamic_model_t::Road_type::ITIME; ++j)
        {
            c_lb[kc] = 0.0;
            c_ub[kc] = 0.0;
            kc.increment();
        }

        // (4.2) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (after time)
        for (size_t j = Dynamic_model_t::Road_type::ITIME+1; j < Dynamic_model_t::NSTATE; ++j)
        {
            c_lb[kc] = 0.0;
            c_ub[kc] = 0.0;
            kc.increment();
        }

        // (4.3) Algebraic constraints: dqa^{i} = 0.0
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
        {
            c_lb[kc] = 0.0;
            c_ub[kc] = 0.0;
            kc.increment();
        }
        
        // (4.4) Extra constraints per point
        auto [c_extra_lb, c_extra_ub] = car.optimal_laptime_extra_constraints_bounds(s[i]); 
        for (size_t j = 0; j < Dynamic_model_t::N_OL_EXTRA_CONSTRAINTS; ++j)
        {
            c_lb[kc] = c_extra_lb[j];
            c_ub[kc] = c_extra_ub[j];
            kc.increment();
        }

        // (4.5) Full mesh control derivatives
        for (size_t j = 0; j < Dynamic_model_t::NCONTROL; ++j)
        {   
            if ( control_variables[j].optimal_control_type == FULL_MESH )
            {
                c_lb[kc] = 0.0;
                c_ub[kc] = 0.0;
                kc.increment();
            }
        }
    }

    // (5) Add the periodic element if track is closed
    if constexpr (isClosed)
    {
        // Equality constraints: --------------- 

        // (5.1) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (before time)
        for (size_t j = 0; j < Dynamic_model_t::Road_type::ITIME; ++j)
        {
            c_lb[kc] = 0.0;
            c_ub[kc] = 0.0;
            kc.increment();
        }

        // (5.2) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (after time)
        for (size_t j = Dynamic_model_t::Road_type::ITIME+1; j < Dynamic_model_t::NSTATE; ++j)
        {
            c_lb[kc] = 0.0;
            c_ub[kc] = 0.0;
            kc.increment();
        }

        // (5.3) Algebraic constraints: dqa^{i} = 0.0
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
        {
            c_lb[kc] = 0.0;
            c_ub[kc] = 0.0;
            kc.increment();
        }
        
        // (5.4) Extra constraints per point
        auto [c_extra_lb, c_extra_ub] = car.optimal_laptime_extra_constraints_bounds(s[0]); 
        for (size_t j = 0; j < Dynamic_model_t::N_OL_EXTRA_CONSTRAINTS; ++j)
        {
            c_lb[kc] = c_extra_lb[j];
            c_ub[kc] = c_extra_ub[j];
            kc.increment();
        }

        // (5.5) Full mesh control variables derivatives
        for (size_t j = 0; j < Dynamic_model_t::NCONTROL; ++j)
        {   
            if ( control_variables[j].optimal_control_type == FULL_MESH )
            {
                c_lb[kc] = 0.0;
                c_ub[kc] = 0.0;
                kc.increment();
            }
        }
    }

    // (6) Check the dimensions of the counter used to fill the vectors
    if ( k != fg.get_n_variables() )
    {
        std::ostringstream s_out;
        s_out << "k(=" << k << ") != fg.get_n_variables()(=" << fg.get_n_variables() << ")" ;
        throw std::runtime_error(s_out.str());
    }
    
    if ( kc != fg.get_n_constraints() )
    {
        std::ostringstream s_out;
        s_out << "kc(=" << kc << ") != fg.get_n_constraints()(=" << fg.get_n_constraints() << ")" ;
        throw std::runtime_error(s_out.str());
    }


    // (7) Warm-start: simply check the dimensions
    if ( warm_start )
    {
        if ( optimization_data.zl.size() != fg.get_n_variables() )
        {
            std::ostringstream s_out;
            s_out << "size of zl should be " << fg.get_n_variables() << " but it is " << optimization_data.zl.size();
            throw std::runtime_error(s_out.str());
        }

        if ( optimization_data.zu.size() != fg.get_n_variables() )
        {
            std::ostringstream s_out;
            s_out << "size of zu should be " << fg.get_n_variables() << " but it is " << optimization_data.zu.size();
            throw std::runtime_error(s_out.str());
        }

        if ( optimization_data.lambda.size() != fg.get_n_constraints() )
        {
            std::ostringstream s_out;
            s_out << "size of lambda should be " << fg.get_n_constraints() << " but it is " << optimization_data.lambda.size();
            throw std::runtime_error(s_out.str());
        }
    } 

    // (8) Run optimization

    // (8.1) Prepare options
    std::string ipoptoptions;
    ipoptoptions += "Integer print_level  ";
    ipoptoptions += std::to_string(options.print_level);
    ipoptoptions += "\n";
    ipoptoptions += "String  sb           yes\n";
    ipoptoptions += "Sparse true forward\n";
    ipoptoptions += "Numeric tol          1e-10\n";
    ipoptoptions += "Numeric constr_viol_tol  1e-10\n";
    ipoptoptions += "Numeric acceptable_tol  1e-8\n";

    // (8.2) Return object
    CppAD::ipopt_cppad_result<std::vector<scalar>> result;

    // (8.3) Solve the problem
    if ( !warm_start )
        CppAD::ipopt_cppad_solve<std::vector<scalar>, FG_derivative<isClosed>>(ipoptoptions, x0, x_lb, x_ub, c_lb, c_ub, fg, result);
    else
        CppAD::ipopt_cppad_solve<std::vector<scalar>, FG_derivative<isClosed>>(ipoptoptions, x0, x_lb, x_ub, c_lb, c_ub, 
            optimization_data.lambda, optimization_data.zl, optimization_data.zu, fg, result);
 

    // (8.4) Check success flag
    success = result.status == CppAD::ipopt_cppad_result<std::vector<scalar>>::success; 
    iter_count = result.iter_count;

    if ( !success && options.throw_if_fail )
    {
        throw std::runtime_error("Optimization did not succeed");
    }

    // (9) Export the solution

    // (9.1) Load the latest solution vector to the fitness function object
    std::vector<CppAD::AD<scalar>> x_final(result.x.size(),0.0),fg_final(fg.get_n_constraints()+1,0.0);
    std::copy(result.x.cbegin(), result.x.cend(), x_final.begin());
    fg(fg_final,x_final);
  
    assert(fg.get_states().size() == n_points);
    assert(fg.get_algebraic_states().size() == n_points);

    // (9.2) Export states
    q = std::vector<std::array<scalar,Dynamic_model_t::NSTATE>>(n_points);
    
    for (size_t i = 0; i < n_points; ++i)
        for (size_t j = 0; j < Dynamic_model_t::NSTATE; ++j)
            q[i][j] = Value(fg.get_state(i)[j]);

    // (9.3) Export algebraic states
    qa = std::vector<std::array<scalar,Dynamic_model_t::NALGEBRAIC>>(n_points);
    
    for (size_t i = 0; i < n_points; ++i)
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
            qa[i][j] = Value(fg.get_algebraic_state(i)[j]);

    // (9.4) Export control variables
    control_variables.clear();

    for (size_t j = 0; j < Dynamic_model_t::NCONTROL; ++j)
    {
        switch (control_variables[j].optimal_control_type) 
        {
         case (CONSTANT):
            std::transform(fg.get_control(j).u.cbegin(), fg.get_control(j).u.cend(), control_variables[j].u.begin(), 
                [](const auto& u_cppad) -> auto { return Value(u_cppad); });
            break;
        
         case (HYPERMESH):
            std::transform(fg.get_control(j).u.cbegin(), fg.get_control(j).u.cend(), control_variables[j].u.begin(), 
                [](const auto& u_cppad) -> auto { return Value(u_cppad); });
            break;

         case (FULL_MESH):
            std::transform(fg.get_control(j).u.cbegin(), fg.get_control(j).u.cend(), control_variables[j].u.begin(), 
                [](const auto& u_cppad) -> auto { return Value(u_cppad); });

            std::transform(fg.get_control(j).dudt.cbegin(), fg.get_control(j).dudt.cend(), control_variables[j].dudt.begin(), 
                [](const auto& u_cppad) -> auto { return Value(u_cppad); });
            break;

         default:
            // do nothing
            break;
        }
    }

    // (9.5) Export time, x, y, and psi
    const scalar& L = car.get_road().track_length();

    x_coord = std::vector<scalar>(fg.get_states().size());
    y_coord = std::vector<scalar>(fg.get_states().size());
    psi = std::vector<scalar>(fg.get_states().size());

    // (9.5.1) Compute the first point
    auto u_i = fg.get_controls().control_array_at_s(car, 0, s.front());
    auto dtimeds_first = fg.get_car()(fg.get_state(0),fg.get_algebraic_state(0),u_i,s[0]).first[Dynamic_model_t::Road_type::ITIME];
    auto dtimeds_prev = dtimeds_first;

    x_coord.front() = Value(fg.get_car().get_road().get_x());
    y_coord.front() = Value(fg.get_car().get_road().get_y());
    psi.front() = Value(fg.get_car().get_road().get_psi());

    // (9.5.2) Compute the rest of the points
    for (size_t i = 1; i < fg.get_states().size(); ++i)
    {
        u_i = fg.get_controls().control_array_at_s(car, i, s[i]);
        const auto dtimeds = fg.get_car()(fg.get_state(i),fg.get_algebraic_state(i),u_i,s[i]).first[Dynamic_model_t::Road_type::ITIME];

        q[i][Dynamic_model_t::Road_type::ITIME] = q[i-1][Dynamic_model_t::Road_type::ITIME] 
            + Value((s[i]-s[i-1])*(options.sigma*dtimeds + (1.0-options.sigma)*dtimeds_prev));
        dtimeds_prev = dtimeds;

        x_coord[i] = Value(fg.get_car().get_road().get_x());
        y_coord[i] = Value(fg.get_car().get_road().get_y());
        psi[i]     = Value(fg.get_car().get_road().get_psi());
    }

    // (9.5.4) Compute the laptime
    laptime = q.back()[Dynamic_model_t::Road_type::ITIME];

    if (isClosed)
        laptime += Value((L-s.back())*((1.0-options.sigma)*dtimeds_prev+options.sigma*dtimeds_first));

    // (9.6) Save optimization data
    optimization_data.zl     = result.zl;
    optimization_data.zu     = result.zu;
    optimization_data.lambda = result.lambda;
}


template<typename Dynamic_model_t>
std::unique_ptr<Xml_document> Optimal_laptime<Dynamic_model_t>::xml() const 
{
    std::unique_ptr<Xml_document> doc_ptr(std::make_unique<Xml_document>());

    const auto [key_name, q_names, qa_names, u_names] = Dynamic_model_t::get_state_and_control_names();

    std::ostringstream s_out;
    s_out.precision(17);

    doc_ptr->create_root_element("optimal_laptime");

    auto root = doc_ptr->get_root_element();
    
    root.add_attribute("n_points",std::to_string(n_points));
    
    if ( is_closed )
        root.add_attribute("type", "closed");
    else
        root.add_attribute("type", "open");

    if ( is_direct )
        root.add_attribute("is_direct","true");
    else
        root.add_attribute("is_direct","false");
    

    // Save arclength
    for (size_t j = 0; j < s.size()-1; ++j)
        s_out << s[j] << ", ";

    s_out << s.back();

    root.add_child("laptime").set_value(std::to_string(laptime));

    root.add_child("arclength").set_value(s_out.str());
    s_out.str(""); s_out.clear();

    // Save state
    for (size_t i = 0; i < Dynamic_model_t::NSTATE; ++i)
    {
        for (size_t j = 0; j < q.size()-1; ++j)
            s_out << q[j][i] << ", " ;

        s_out << q.back()[i];

        root.add_child(q_names[i]).set_value(s_out.str());
        s_out.str(""); s_out.clear();
    }

    // Save algebraic state
    for (size_t i = 0; i < Dynamic_model_t::NALGEBRAIC; ++i)
    {
        for (size_t j = 0; j < qa.size()-1; ++j)
            s_out << qa[j][i] << ", " ;

        s_out << qa.back()[i];

        root.add_child(qa_names[i]).set_value(s_out.str());
        s_out.str(""); s_out.clear();
    }

    // Save controls
    auto node_control_variables = root.add_child("control_variables");
    for (size_t i = 0; i < Dynamic_model_t::NCONTROL; ++i)
    {
        auto node_variable = node_control_variables.add_child(u_names[i]);
        if ( control_variables[i].optimal_control_type == DONT_OPTIMIZE )
        {
            node_variable.add_attribute("type","dont optimize");
            node_variable.add_child("values");
        }
        else if ( control_variables[i].optimal_control_type == CONSTANT )
        {   
            node_variable.add_attribute("type","constant");
            s_out << control_variables[i].u.front();
            node_variable.add_child("values").set_value(s_out.str());
            s_out.str(""); s_out.clear();
        }
        else if ( control_variables[i].optimal_control_type == HYPERMESH )
        {
            node_variable.add_attribute("type","hypermesh");
            for (size_t j = 0; j < control_variables[i].u.size()-1; ++j)
                s_out << control_variables[i].u[j] << ", " ;

            s_out << control_variables[i].u.back();

            node_variable.add_child("values").set_value(s_out.str());
            s_out.str(""); s_out.clear();

            for (size_t j = 0; j < control_variables[i].s_hypermesh.size()-1; ++j)
                s_out << control_variables[i].s_hypermesh[j] << ", " ;

            s_out << control_variables[i].s_hypermesh.back();

            node_variable.add_child("hypermesh").set_value(s_out.str());
            s_out.str(""); s_out.clear();
        }
        else if ( control_variables[i].optimal_control_type == FULL_MESH )
        {
            node_variable.add_attribute("type", "full-mesh");
            for (size_t j = 0; j < control_variables[i].u.size()-1; ++j)
                s_out << control_variables[i].u[j] << ", " ;

            s_out << control_variables[i].u.back();

            node_variable.add_child("values").set_value(s_out.str());
            s_out.str(""); s_out.clear();

            // Save controls derivatives if not direct
            if ( !is_direct )
            {
                for (size_t j = 0; j < control_variables[i].dudt.size()-1; ++j)
                    s_out << control_variables[i].dudt[j] << ", " ;
        
                s_out << control_variables[i].dudt.back();

                node_variable.add_child("derivatives").set_value(s_out.str());
                s_out.str(""); s_out.clear();
            }
        }
    }

    // Save x, y, and psi if they are not contained
    for (size_t j = 0; j < q.size()-1; ++j)
        s_out << x_coord[j] << ", " ;

    s_out << x_coord.back();
    root.add_child("x").set_value(s_out.str());
    s_out.str(""); s_out.clear();

    for (size_t j = 0; j < q.size()-1; ++j)
        s_out << y_coord[j] << ", " ;

    s_out << y_coord.back();
    root.add_child("y").set_value(s_out.str());
    s_out.str(""); s_out.clear();

    for (size_t j = 0; j < q.size()-1; ++j)
        s_out << psi[j] << ", " ;

    s_out << psi.back();
    root.add_child("psi").set_value(s_out.str());
    s_out.str(""); s_out.clear();

    // Save optimization data
    auto opt_data = root.add_child("optimization_data");
    size_t n_vars(0);
    size_t n_cons(0);

    if ( is_direct )
    {
        n_vars = n_variables_per_point<true>(control_variables);
        n_cons = n_constraints_per_element<true>(control_variables);
    }
    else
    {
        n_vars = n_variables_per_point<false>(control_variables);
        n_cons = n_constraints_per_element<false>(control_variables);
    }

    opt_data.add_attribute("n_variables_per_point", std::to_string(n_vars));
    opt_data.add_attribute("n_constraints_per_element", std::to_string(n_cons));

    for (size_t j = 0; j < optimization_data.zl.size()-1; ++j)
        s_out << optimization_data.zl[j] << ", " ;

    s_out << optimization_data.zl.back();

    opt_data.add_child("zl").set_value(s_out.str());
    s_out.str(""); s_out.clear();

    for (size_t j = 0; j < optimization_data.zu.size()-1; ++j)
        s_out << optimization_data.zu[j] << ", " ;

    s_out << optimization_data.zu.back();

    opt_data.add_child("zu").set_value(s_out.str());
    s_out.str(""); s_out.clear();

    for (size_t j = 0; j < optimization_data.lambda.size()-1; ++j)
        s_out << optimization_data.lambda[j] << ", " ;

    s_out << optimization_data.lambda.back();

    opt_data.add_child("lambda").set_value(s_out.str());
    s_out.str(""); s_out.clear();

    return doc_ptr;
}


template<typename Dynamic_model_t>
template<bool isClosed>
inline void Optimal_laptime<Dynamic_model_t>::FG_direct<isClosed>::operator()(FG_direct<isClosed>::ADvector& fg, const Optimal_laptime<Dynamic_model_t>::FG_direct<isClosed>::ADvector& x)
{
    // (1) Useful aliases
    auto& _n_points          = FG::_n_points;
    auto& _car               = FG::_car;
    auto& _s                 = FG::_s;
    auto& _q0                = FG::_q0;
    auto& _qa0               = FG::_qa0;
    auto& _u0                = FG::_u0;
    auto& _q                 = FG::_q;
    auto& _qa                = FG::_qa;
    auto& _control_variables = FG::_control_variables;
    auto& _dqdt              = FG::_dqdt;
    auto& _dqa               = FG::_dqa ;
    auto& _sigma             = FG::_sigma;

    // (2) Check dimensions
    assert(x.size() == FG::_n_variables);
    assert(fg.size() == (1 + FG::_n_constraints));

    // (3) Load the state and control vectors
    size_t k = 0;

    // (3.1) For open simulations, load the stored initial condition
    if constexpr (!isClosed)
    {
        // (3.1.1) State
        for (size_t j = 0; j < Dynamic_model_t::NSTATE; ++j)
            _q[0][j] = _q0[j];

        // (3.1.2) Algebraic state
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
            _qa[0][j] = _qa0[j];
    
        // (3.1.3) Control for full mesh calculations
        for (size_t j = 0; j < Dynamic_model_t::NCONTROL; ++j)
        {
            if ( _control_variables[j].optimal_control_type == FULL_MESH )
                _control_variables[j].u.front() = _u0[j];
        }
    }

    // (3.2) Load the rest of the points
    constexpr const size_t offset = (isClosed ? 0 : 1);

    for (size_t i = offset; i < _n_points; ++i)
    {
        // (3.2.1) Load state (before time)
        for (size_t j = 0; j < Dynamic_model_t::Road_type::ITIME; ++j)
            _q[i][j] = x[k++];

        // (3.2.2) Load state (after time)
        for (size_t j = Dynamic_model_t::Road_type::ITIME+1; j < Dynamic_model_t::NSTATE; ++j)
            _q[i][j] = x[k++];

        // (3.2.3) Load algebraic state
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
            _qa[i][j] = x[k++];
            
        // (3.2.4) Load full mesh controls
        for (auto& control_variable : _control_variables)
        {
            if ( control_variable.optimal_control_type == FULL_MESH )
                control_variable.u[i] = x[k++];
        }
    }

    // (3.3) Load the rest of the controls
    for (auto& control_variable : _control_variables)
    {
        switch (control_variable.optimal_control_type)
        {
         case (CONSTANT):
            for (size_t i = 0; i < control_variable.u.size(); ++i)
                control_variable.u[i] = x[k++];
            break;
         case (HYPERMESH):
            for (size_t i = 0; i < control_variable.u.size(); ++i)
                control_variable.u[i] = x[k++];
            break;
         default:
            break;
        }
    }

    // (4) Check that all variables in x were used
    assert(k == FG::_n_variables);

    // (5) Write fitness function and constraints
    fg[0] = 0.0;
    
    std::tie(_dqdt[0], _dqa[0]) = _car(_q[0],_qa[0],_control_variables.control_array_at_s(_car, 0,_s[0]),_s[0]);
    k = 1;  // Reset the counter
    for (size_t i = 1; i < _n_points; ++i)
    {
        std::tie(_dqdt[i],_dqa[i]) = _car(_q[i],_qa[i],_control_variables.control_array_at_s(_car, i, _s[i]),_s[i]);

        // (5.1) Fitness function: integral of time
        //fg[0] += (_s[i]-_s[i-1])*((1.0-_sigma)*_dqdt[i-1][Dynamic_model_t::Road_type::ITIME] + _sigma*_dqdt[i][Dynamic_model_t::Road_type::ITIME]);
        fg[0] += (_s[i]-_s[i-1])*(_q[i][0]-205.777)*(_q[i][0]-205.777);
            

        // Equality constraints: --------------- 

        // (5.2) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (before time)
        for (size_t j = 0; j < Dynamic_model_t::Road_type::ITIME; ++j)
            fg[k++] = _q[i][j] - _q[i-1][j] - (_s[i]-_s[i-1])*((1.0-_sigma)*_dqdt[i-1][j] + _sigma*_dqdt[i][j]);

        // (5.3) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (after time)
        for (size_t j = Dynamic_model_t::Road_type::ITIME+1; j < Dynamic_model_t::NSTATE; ++j)
            fg[k++] = _q[i][j] - _q[i-1][j] - (_s[i]-_s[i-1])*((1.0-_sigma)*_dqdt[i-1][j] + _sigma*_dqdt[i][j]);

        // (5.4) algebraic constraints: dqa^{i} = 0.0
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
            fg[k++] = _dqa[i][j];
        
        // (5.5) Problem dependent extra constraints
        auto c_extra = _car.optimal_laptime_extra_constraints();

        for (size_t j = 0; j < Dynamic_model_t::N_OL_EXTRA_CONSTRAINTS; ++j)
            fg[k++] = c_extra[j];
    }

    // (5.6) Add a penalisation to the controls
    for (const auto& control_variable : _control_variables)
    {
        if (control_variable.optimal_control_type == FULL_MESH)
        {
            for (size_t i = 1; i < _n_points; ++i)
            {
                const auto derivative = (control_variable.u[i]-control_variable.u[i-1])/(_s[i]-_s[i-1]);
                fg[0] += control_variable.dissipation*(derivative*derivative)*(_s[i]-_s[i-1]);
            }
        }
    }

    // (5.7) Add the periodic element if track is closed
    const scalar& L = _car.get_road().track_length();
    if constexpr (isClosed)
    {
        // (5.7.1) Fitness function: integral of time
        fg[0] += (L-_s.back())*(_sigma*_dqdt.front()[Dynamic_model_t::Road_type::ITIME] + (1.0-_sigma)*_dqdt.back()[Dynamic_model_t::Road_type::ITIME]);

        // Equality constraints: 

        // (5.7.2) q^{0} = q^{n-1} + 0.5.ds.[dqdt^{0} + dqdt^{n-1}] (before time)
        for (size_t j = 0; j < Dynamic_model_t::Road_type::ITIME; ++j)
            fg[k++] = _q.front()[j] - _q.back()[j] - (L-_s.back())*((1.0-_sigma)*_dqdt.back()[j] + _sigma*_dqdt.front()[j]);

        // (5.7.3) q^{0} = q^{n-1} + 0.5.ds.[dqdt^{0} + dqdt^{n-1}] (after time)
        for (size_t j = Dynamic_model_t::Road_type::ITIME+1; j < Dynamic_model_t::NSTATE; ++j)
            fg[k++] = _q.front()[j] - _q.back()[j] - (L-_s.back())*((1.0-_sigma)*_dqdt.back()[j] + _sigma*_dqdt.front()[j]);

        // (5.7.4) algebraic constraints: dqa^{0} = 0.0
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
            fg[k++] = _dqa.front()[j];

        // (5.7.5) Inequality constraints: -0.11 < kappa < 0.11, -0.11 < lambda < 0.11
        _car(_q[0],_qa[0],_control_variables.control_array_at_s(_car,0,_s[0]),_s[0]);
        const auto c_extra = _car.optimal_laptime_extra_constraints();
    
        for (size_t j = 0; j < Dynamic_model_t::N_OL_EXTRA_CONSTRAINTS; ++j)
            fg[k++] = c_extra[j];

        // (5.7.6) Add a penalisation to the controls
        for (const auto& control_variable : _control_variables)
        {
            if (control_variable.optimal_control_type == FULL_MESH)
            {
                const auto derivative = (control_variable.u.front()-control_variable.u.back())/(L-_s.back());
                fg[0] += control_variable.dissipation*(derivative*derivative)*(L-_s.back());
            }
        }
    }

    assert(k == FG::_n_constraints+1);
}


template<typename Dynamic_model_t>
template<bool isClosed>
inline void Optimal_laptime<Dynamic_model_t>::FG_derivative<isClosed>::operator()(FG_derivative<isClosed>::ADvector& fg, const Optimal_laptime<Dynamic_model_t>::FG_derivative<isClosed>::ADvector& x)
{
    // (1) Useful aliases
    auto& _n_points          = FG::_n_points;
    auto& _car               = FG::_car;
    auto& _s                 = FG::_s;
    auto& _q0                = FG::_q0;
    auto& _qa0               = FG::_qa0;
    auto& _u0                = FG::_u0;
    auto& _q                 = FG::_q;
    auto& _qa                = FG::_qa;
    auto& _control_variables = FG::_control_variables;
    auto& _dqdt              = FG::_dqdt;
    auto& _dqa               = FG::_dqa ;
    auto& _sigma             = FG::_sigma;

    // (2) Check dimensions
    assert(x.size() == FG::_n_variables);
    assert(fg.size() == (1 + FG::_n_constraints));

    // (3) Load the state and control vectors
    size_t k = 0;

    // (3.1) For open simulations, load the stored initial condition
    if constexpr (!isClosed)
    {
        // (3.1.1) State
        for (size_t j = 0; j < Dynamic_model_t::NSTATE; ++j)
            _q[0][j] = _q0[j];

        // (3.1.2) Algebraic state
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
            _qa[0][j] = _qa0[j];
    
        // (3.1.3) Control for full mesh calculations
        for (size_t j = 0; j < Dynamic_model_t::NCONTROL; ++j)
        {
            if ( _control_variables[j].optimal_control_type == FULL_MESH )
            {
                _control_variables[j].u.front()    = _u0[j];
                _control_variables[j].dudt.front() = _dudt0[j];
            }
        }
    }

    // (3.2) Load the rest of the points
    constexpr const size_t offset = (isClosed ? 0 : 1);

    for (size_t i = offset; i < _n_points; ++i)
    {
        // (3.2.1) Load state (before time)
        for (size_t j = 0; j < Dynamic_model_t::Road_type::ITIME; ++j)
            _q[i][j] = x[k++];

        // (3.2.2) Load state (after time)
        for (size_t j = Dynamic_model_t::Road_type::ITIME+1; j < Dynamic_model_t::NSTATE; ++j)
            _q[i][j] = x[k++];

        // (3.2.3) Load algebraic state
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
            _qa[i][j] = x[k++];
            
        // (3.2.4) Load full mesh controls
        for (auto& control_variable : _control_variables)
        {
            if ( control_variable.optimal_control_type == FULL_MESH )
                control_variable.u[i] = x[k++];
        }

        // (3.2.4) Load full mesh control derivatives
        for (auto& control_variable : _control_variables)
        {
            if ( control_variable.optimal_control_type == FULL_MESH )
                control_variable.dudt[i] = x[k++];
        }
    }

    // (3.3) Load the rest of the controls
    for (auto& control_variable : _control_variables)
    {
        switch (control_variable.optimal_control_type)
        {
         case (CONSTANT):
            for (size_t i = 0; i < control_variable.u.size(); ++i)
                control_variable.u[i] = x[k++];
            break;
         case (HYPERMESH):
            for (size_t i = 0; i < control_variable.u.size(); ++i)
                control_variable.u[i] = x[k++];
            break;
         default:
            break;
        }
    }


    // (4) Check that all variables in x were used
    assert(k == FG::_n_variables);

    // (5) Write fitness function and constraints
    fg[0] = 0.0;

    std::tie(_dqdt[0], _dqa[0]) = _car(_q[0],_qa[0],_control_variables.control_array_at_s(_car,0,_s[0]),_s[0]);
    k = 1;  // Reset the counter
    for (size_t i = 1; i < _n_points; ++i)
    {
        std::tie(_dqdt[i],_dqa[i]) = _car(_q[i],_qa[i],_control_variables.control_array_at_s(_car,i,_s[i]),_s[i]);

        // (5.1) Fitness function: 

        // (5.1.1) Integral of time
        fg[0] += (_s[i]-_s[i-1])*((1.0-_sigma)*_dqdt[i-1][Dynamic_model_t::Road_type::ITIME] + _sigma*_dqdt[i][Dynamic_model_t::Road_type::ITIME]);

        // (5.1.2) Penalisation of controls
        for (const auto& control_variable : _control_variables)
        {
            if ( control_variable.optimal_control_type == FULL_MESH )
            {
                const auto& dudt = control_variable.dudt;
                const auto& dissipation = control_variable.dissipation;
                fg[0] += 0.5*dissipation*(dudt[i-i]*dudt[i-i]+dudt[i]*dudt[i])*(_s[i]-_s[i-1]);
            }
        }

        // Equality constraints: 

        // (5.2) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (before time)
        for (size_t j = 0; j < Dynamic_model_t::Road_type::ITIME; ++j)
            fg[k++] = _q[i][j] - _q[i-1][j] - (_s[i]-_s[i-1])*((1.0-_sigma)*_dqdt[i-1][j] + _sigma*_dqdt[i][j]);

        // (5.3) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (before time)
        for (size_t j = Dynamic_model_t::Road_type::ITIME+1; j < Dynamic_model_t::NSTATE; ++j)
            fg[k++] = _q[i][j] - _q[i-1][j] - (_s[i]-_s[i-1])*((1.0-_sigma)*_dqdt[i-1][j] + _sigma*_dqdt[i][j]);

        // (5.4) algebraic constraints: dqa^{i} = 0.0
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
            fg[k++] = _dqa[i][j];

        // (5.5) Problem dependent extra constraints
        auto c_extra = _car.optimal_laptime_extra_constraints();

        for (size_t j = 0; j < Dynamic_model_t::N_OL_EXTRA_CONSTRAINTS; ++j)
            fg[k++] = c_extra[j];

        // (5.6) Time derivative of full mesh controls
        for (const auto& control_variable : _control_variables)
        {
            if (control_variable.optimal_control_type == FULL_MESH)
            {
                const auto& u = control_variable.u;
                const auto& dudt = control_variable.dudt;
                fg[k++] = u[i] - u[i-1] - (_s[i]-_s[i-1])
                    *((1.0-_sigma)*dudt[i-1]*_dqdt[i-1][Dynamic_model_t::Road_type::ITIME]+_sigma*dudt[i]*_dqdt[i][Dynamic_model_t::Road_type::ITIME]);
            }
        }
    }

    // (5.7) Add the periodic element if track is closed
    const scalar& L = _car.get_road().track_length();
    if constexpr (isClosed)
    {
        // (5.7.1) Fitness function: 

        // (5.7.1.1) Integral of time
        fg[0] += (L-_s.back())*((1.0-_sigma)*_dqdt.back()[Dynamic_model_t::Road_type::ITIME] + _sigma*_dqdt.front()[Dynamic_model_t::Road_type::ITIME]);

        // (5.7.1.2) Penalisation to the controls
        for (const auto& control_variable : _control_variables)
        {
            if ( control_variable.optimal_control_type == FULL_MESH )
            {
                const auto& dudt = control_variable.dudt;
                const auto& dissipation = control_variable.dissipation;
                fg[0] += 0.5*dissipation*(dudt.back()*dudt.back() + dudt.front()*dudt.front())*(L-_s.back());
            }
        }

        // Equality constraints: 

        // (5.7.2) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (before time)
        for (size_t j = 0; j < Dynamic_model_t::Road_type::ITIME; ++j)
            fg[k++] = _q.front()[j] - _q.back()[j] - (L-_s.back())*((1.0-_sigma)*_dqdt.back()[j] + _sigma*_dqdt.front()[j]);

        // (5.7.3) q^{i} = q^{i-1} + 0.5.ds.[dqdt^{i} + dqdt^{i-1}] (before time)
        for (size_t j = Dynamic_model_t::Road_type::ITIME+1; j < Dynamic_model_t::NSTATE; ++j)
            fg[k++] = _q.front()[j] - _q.back()[j] - (L-_s.back())*((1.0-_sigma)*_dqdt.back()[j] + _sigma*_dqdt.front()[j]);

        // (5.7.3) algebraic constraints: dqa^{0} = 0
        for (size_t j = 0; j < Dynamic_model_t::NALGEBRAIC; ++j)
            fg[k++] = _dqa.front()[j];

        // (5.7.4) Inequality constraints: -0.11 < kappa < 0.11, -0.11 < lambda < 0.11
        _car(_q[0],_qa[0],_control_variables.control_array_at_s(_car,0,_s[0]),_s[0]);
        auto c_extra = _car.optimal_laptime_extra_constraints();

        for (size_t j = 0; j < Dynamic_model_t::N_OL_EXTRA_CONSTRAINTS; ++j)
            fg[k++] = c_extra[j];

        // (5.7.5) Time derivative of the control variables
        for (const auto& control_variable : _control_variables)
        {
            if (control_variable.optimal_control_type == FULL_MESH)
            {
                const auto& u = control_variable.u;
                const auto& dudt = control_variable.dudt;
                fg[k++] = u.front() - u.back() - (L-_s.back())
                    *((1.0-_sigma)*dudt.back()*_dqdt.back()[Dynamic_model_t::Road_type::ITIME]+_sigma*dudt.front()*_dqdt.front()[Dynamic_model_t::Road_type::ITIME]);
            }
        }
    }

    assert(k == FG::_n_constraints+1);
}

//
// Implementation of the helper class Control_variables ----------------------------------------------------:-
//

template<typename Dynamic_model_t>
template<typename T>
void Optimal_laptime<Dynamic_model_t>::Control_variables<T>::check_inputs()
{
    for (auto& control_variable : *this)
    {
        switch(control_variable.optimal_control_type)
        {
         case(DONT_OPTIMIZE):
            // Check that s_hypermesh is empty, and that u only contains one value 
            if ( control_variable.s_hypermesh.size() > 0 ) 
                throw std::runtime_error("[ERROR] Control_variables::check_inputs() -> "
                    "In \"don't optimize\" mode, s_hypermesh should be empty");

            if ( control_variable.u.size() > 0 )
                throw std::runtime_error("[ERROR] Control_variables::check_inputs() -> "
                    "In \"don't optimize\" mode, u should be empty. Its size is " + std::to_string(control_variable.u.size()) );

            if ( control_variable.dudt.size() != 0 )
                throw std::runtime_error("[ERROR] Control_variables::check_inputs() -> "
                    "In \"don't optimize\" mode, dudt should be empty");

            control_variable.dissipation = 0.0;

            break;

         case(CONSTANT):
            // Check that s_hypermesh is empty, and that u only contains one value 
            if ( control_variable.s_hypermesh.size() > 0 ) 
                throw std::runtime_error("[ERROR] Control_variables::check_inputs() -> "
                    "In \"constant optimization\" mode, s_hypermesh should be empty");

            if ( control_variable.u.size() != 1 )
                throw std::runtime_error("[ERROR] Control_variables::check_inputs() -> "
                    "In \"constant optimization\" mode, u should contain only one value");

            if ( control_variable.dudt.size() != 0 )
                throw std::runtime_error("[ERROR] Control_variables::check_inputs() -> "
                    "In \"constant optimization\" mode, dudt should be empty");

            control_variable.dissipation = 0.0;

            break;
         case(HYPERMESH):
            // Check that the size of s_hypermesh equals the size of u
            if ( control_variable.s_hypermesh.size() != control_variable.u.size() )
                throw std::runtime_error("[ERROR] Control_variables::check_inputs() -> "
                    "In \"hypermesh optimization\" mode, size(s_hypermesh) should be equal to size(u)");

            if ( control_variable.dudt.size() != 0 )
                throw std::runtime_error("[ERROR] Control_variables::check_inputs() -> "
                    "In \"hypermesh optimization\" mode, dudt should be empty");

            control_variable.dissipation = 0.0;

            break;
         case(FULL_MESH): 
            // Check that s_hypermesh is empty, the size of u will be checked later 
            if ( control_variable.s_hypermesh.size() > 0 ) 
                throw std::runtime_error("[ERROR] Control_variables::check_inputs() -> "
                    "In \"full-mesh optimize\" mode, s_hypermesh should be empty");

            if ( (control_variable.u.size() != control_variable.dudt.size()) && (control_variable.dudt.size() != 0) )
                throw std::runtime_error("[ERROR] Control_variables::check_inputs() -> "
                    "In \"full-mesh optimize\" mode, dudt should either have the same size as u, or be empty");

            if ( control_variable.dissipation < 0.0 )
                throw std::runtime_error("[ERROR] Control_variables::check_inputs() -> dissipation must be non-negative");

            break;
         default:
            throw std::runtime_error("[ERROR] Control_variables::check_inputs() -> optimization mode not recognized");
    
            break;
        }
    }
}


template<typename Dynamic_model_t>
template<typename T>
void Optimal_laptime<Dynamic_model_t>::Control_variables<T>::compute_statistics()
{
    // (1) Reset values
    number_of_constant_optimizations        = 0;
    number_of_hypermesh_optimization_points = 0;
    number_of_full_optimizations            = 0;

    for (const auto& control_variable : *this)
    {
        switch(control_variable.optimal_control_type)
        {
         case(DONT_OPTIMIZE):
            // Do nothing
            break;
         case(CONSTANT):
            ++number_of_constant_optimizations;
            break;
         case(HYPERMESH):
            number_of_hypermesh_optimization_points += control_variable.s_hypermesh.size();
            break;
         case(FULL_MESH): 
            ++number_of_full_optimizations;
            break;
         default:
            throw std::runtime_error("[ERROR] Control_variables::check_inputs() -> optimization mode not recognized");
    
            break;
        }
    }
}


template<typename Dynamic_model_t>
template<typename T>
template<typename U>
std::enable_if_t<std::is_same_v<U,scalar>, typename Optimal_laptime<Dynamic_model_t>::template Control_variables<CppAD::AD<scalar>>> 
    Optimal_laptime<Dynamic_model_t>::Control_variables<T>::to_CppAD() const
{
    // (1) Create return value
    Control_variables<CppAD::AD<scalar>> output;

    // (2) Copy the std::array storage
    std::transform(base_type::cbegin(), base_type::cend(), output.begin(), [](const auto& scalar_control) -> auto { return scalar_control.to_CppAD(); });

    // (3) Copy the statistics
    output.number_of_constant_optimizations        = number_of_constant_optimizations;
    output.number_of_hypermesh_optimization_points = number_of_hypermesh_optimization_points;
    output.number_of_full_optimizations            = number_of_full_optimizations;

    // (4) Return
    return output;
}


template<typename Dynamic_model_t>
template<typename T>
template<typename U>
std::enable_if_t<std::is_same_v<U,scalar>, typename Optimal_laptime<Dynamic_model_t>::template Control_variable<CppAD::AD<scalar>>> 
    Optimal_laptime<Dynamic_model_t>::Control_variable<T>::to_CppAD() const
{
    // (1) Define the new variable of type CppAD::AD
    Control_variable<CppAD::AD<scalar>> output;

    // (2) Copy optimal control type
    output.optimal_control_type = optimal_control_type; 

    // (3) Copy hypermesh
    output.s_hypermesh          = s_hypermesh; 

    // (4) Copy control variable values
    output.u = std::vector<CppAD::AD<scalar>>(u.size());
    std::copy(u.cbegin(), u.cend(), output.u.begin());

    output.dudt = std::vector<CppAD::AD<scalar>>(dudt.size());
    std::copy(dudt.cbegin(), dudt.cend(), output.dudt.begin());

    // (5) Copy dissipation
    output.dissipation = dissipation;

    // (6) Return
    return output;
}
