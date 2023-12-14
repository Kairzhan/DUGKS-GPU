#ifndef CONFIG_H_
#define CONFIG_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>

using namespace std;

//***************************************************************************
// Simple class to read configuration from the text file. Note
// that input filename is hardcoded as config.txt
//***************************************************************************
class Config {
public:
    Config()
    {
        std::vector<std::string> lines;
        std::string temp;
        std::ifstream file("config.txt");

        while (std::getline(file, temp, '\n')) {
            lines.push_back(temp);
        }

        for (const auto& s : lines) {
            std::vector<string> v;

            split(s, '=', v);

            if (v.size()>1)
                params.insert({v[0], v[1]});
        }
    }

    //===========================================================================
    // Returns string value of the input parameter "k"
    //===========================================================================
    string get(string k)
    {
        std::cout << "reading..." << k << ": " << params[k] << std::endl;
        try {
            return params[k];
        } catch (...) {
            std::cout << "key not found: " << k << std::endl;
            std:exit(1);
        }
    }

private:
    
    std::map<std::string, std::string> params;

    // orig. source:https://www.oreilly.com/library/view/c-cookbook/0596007612/ch04s07.html
    void split(const string& s, char c,
               vector<string>& v)
    {
        string::size_type i = 0;
        string::size_type j = s.find(c);

        while (j != string::npos) {
            v.push_back(s.substr(i, j-i));
            i = ++j;
            j = s.find(c, j);

            if (j == string::npos)
                v.push_back(s.substr(i, s.length()));
        }
    }

};

#endif
