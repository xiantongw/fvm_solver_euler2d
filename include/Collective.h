//
// Created by Xiantong Wang on 2019-04-10.
//

#ifdef _MSC_VER 
//not #if defined(_WIN32) || defined(_WIN64) because we have strncasecmp in mingw
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

#ifndef COLLECTIVE_H
#define COLLECTIVE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

#include "../include/Param.h"

Param ReadParamIn(std::string param_filename);

#endif
