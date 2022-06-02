#ifndef TR7_STRING_H
#define TR7_STRING_H
#pragma once


#include <string>
#include <cctype>
#include <algorithm>

#include "tr7/tr7defs.h"


//
// Defines general utilities to deal with std::strings.
//


namespace tr7
{
    inline pure_function void str2upper_in_place(std::string& str)
    {
        //
        // Converts the input std::string to upper case,
        // in-place.
        //

        std::for_each(str.begin(), str.end(),
                      [](char &c) { c = static_cast<char>(std::toupper(static_cast<unsigned char>(c))); });
    }

    inline pure_function std::string str2upper(std::string str)
    {
        //
        // Converts the input std::string to upper case
        // and returns the result.
        //

        str2upper_in_place(str);
        return str;
    }


    inline pure_function void str2lower_in_place(std::string& str)
    {
        //
        // Converts the input std::string to lower case,
        // in-place.
        //

        std::for_each(str.begin(), str.end(),
            [](char &c) { c = static_cast<char>(std::tolower(static_cast<unsigned char>(c))); });
    }

    inline pure_function std::string str2lower(std::string str)
    {
        //
        // Converts the input std::string to lower case
        // and returns the result.
        //

        str2lower_in_place(str);
        return str;
    }


    inline pure_function void strrep_in_place(std::string& str, const std::string& from, const std::string& to)
    {
        //
        // Replaces all occurrences of "from" with "to" in the
        // input std::string "str", in place.
        //

        if (from.empty()) {
            return;
        }
        std::size_t start_pos{ 0u };
        while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
            str.replace(start_pos, from.length(), to);
            start_pos += to.length(); // in case 'to' contains 'from', like replacing 'x' with 'yx'
        }

    }

    inline pure_function std::string strrep(std::string str, const std::string& from, const std::string& to)
    {
        //
        // Replaces all occurrences of "from" with "to" in the input
        // std::string and returns the result.
        //

        strrep_in_place(str, from, to);
        return str;
    }


    inline pure_function void ltrim_in_place(std::string &s)
    {
        //
        // Left trims an std::string, in-place.
        //

        s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            [](unsigned char ch)
        {
            return !std::isspace(ch);
        }));

    }

    inline pure_function void rtrim_in_place(std::string &s)
    {
        //
        // Right trims an std::string, in-place.
        //

        s.erase(std::find_if(s.rbegin(), s.rend(),
            [](unsigned char ch)
        {
            return !std::isspace(ch);
        }).base(), s.end());

    }

    inline pure_function void trim_in_place(std::string &s)
    {
        //
        // Trims an std::string, in-place.
        //

        ltrim_in_place(s);
        rtrim_in_place(s);
    }


    inline pure_function std::string ltrim(std::string s)
    {
        //
        // Left-trims an std::string, returns the trimmed string.
        //

        ltrim_in_place(s);
        return s;
    }

    inline pure_function std::string rtrim(std::string s)
    {
        //
        // Right-trims an std::string, returns the trimmed string.
        //

        rtrim_in_place(s);
        return s;
    }

    inline pure_function std::string trim(std::string s)
    {
        //
        // Trims an std::string, returns the trimmed string.
        //

        trim_in_place(s);
        return s;
    }


}


#endif
