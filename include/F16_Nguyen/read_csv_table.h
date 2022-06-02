#ifndef F16_NGUYEN_READ_CSV_TABLE_H
#define F16_NGUYEN_READ_CSV_TABLE_H
#pragma once


#include <fstream>
#include <sstream>
#include <iterator>
#include <tuple>
#include <vector>

#include "tr7/tr7string.h"


//
// Defines the "read_csv_table" function, that reads
// the dataset files we have for the F16 model.
//


namespace F16_Nguyen
{
    template<class CharT, class Traits, class Allocator>
    inline std::basic_istream<CharT, Traits>& safe_getline(std::basic_istream<CharT, Traits>& input,
                                                           std::basic_string<CharT, Traits, Allocator> &str,
                                                           const std::vector<std::string> &comments);


    template<typename T = double>
    inline std::tuple<std::vector<T>, std::size_t, std::size_t> read_csv_table(const std::string &csv_filename,
                                                                               char separator = ',',
                                                                               const std::vector<std::string> &comment = { "#" })
    {
        //
        // Returns an std::tuple holding "[table_colmajor, num_rows, num_columns]",
        // read from a "*.csv" source.
        //

        std::ifstream file(csv_filename);
        if (!file.is_open()) {
            TR7_THROW_RUNTIME_ERROR(std::string{ "F16_Nguyen::read_csv_table: cannot open file" } +
                                                 " \"" + csv_filename + "\".")

        }

        std::size_t num_columns{ 0u };
        std::size_t num_rows{ 0u };
        std::vector<T> table_rowmaj;
        for (std::string line; safe_getline(file, line, comment) && !file.eof(); ) {
            std::istringstream ss(line);
            std::size_t line_columns{ 0u };
            for (T t; ss >> t; ) {
                ++line_columns;
                table_rowmaj.push_back(t);
                if (ss.peek() == separator) {
                    ss.ignore();
                }

            }

            if (num_columns == 0u) {
                num_columns = line_columns;

            }
            else if (line_columns != num_columns) {
                TR7_THROW_RUNTIME_ERROR(std::string{ "F16_Nguyen::read_csv_table: inconsistent \"num_columns\"" } +
                                                     " while reading file \"" + csv_filename + "\".")

            }

            ++num_rows;

        }

        // transpose the table to get it in column-major storage
        std::vector<T> table_colmaj(num_rows * num_columns);
        for (auto col = 0u; col < num_columns; ++col) {
            for (auto row = 0u; row < num_rows; ++row) {
                table_colmaj[row + num_rows * col] = table_rowmaj[col + num_columns * row];
            }

        }

        // return the result
        return std::make_tuple(table_colmaj, num_rows, num_columns);

    }


    template<class CharT, class Traits, class Allocator>
    inline std::basic_istream<CharT, Traits>& safe_getline(std::basic_istream<CharT, Traits>& input,
                                                           std::basic_string<CharT, Traits, Allocator> &str,
                                                           const std::vector<std::string> &comments)
    {
        //
        // Works like std::getline, but ignoring empty or commented lines.
        //

        const auto is_empty_or_comment = [&comments](const std::string &trimmed_str)
        {
            if (!trimmed_str.empty()) {
                // does the string start (from the first non-whitespace character) with a comment?
                constexpr auto strcmp_success = 0;
                for (const auto &c : comments) {
                    if (trimmed_str.compare(0, std::min(trimmed_str.size(), c.size()), c) == strcmp_success) {
                        return true;
                    }
                }

                return false;

            }
            else {
                return true;
            }

        };

        do {
            std::getline(input, str);
            tr7::trim_in_place(str);

        } while (is_empty_or_comment(str) && !input.eof() && input);

        if ((input.eof() || !input) && is_empty_or_comment(str)) {
            // for files that end with a comment...
            str.clear();

        }

        return input;

    }


}


#endif
