/** This header contains the definitions of the exceptions that will be used for
  * error handling in this program.
  */

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <stdexcept>

namespace Amber {

class AmberParmError : public std::runtime_error {
    public:
        AmberParmError(std::string const& s) :
            std::runtime_error(s) {}
        AmberParmError(const char* s) :
            std::runtime_error(std::string(s)) {}
};

class AmberCrdError : public std::runtime_error {
    public:
        AmberCrdError(std::string const& s) :
            std::runtime_error(s) {}
        AmberCrdError(const char* s) :
            std::runtime_error(std::string(s)) {}
};

class NotNetcdf : public std::runtime_error {
    public:
        NotNetcdf(std::string const& s) :
            std::runtime_error(s) {}
        NotNetcdf(const char* s) :
            std::runtime_error(std::string(s)) {}
};

class UnitCellError : public std::runtime_error {
    public:
        UnitCellError(std::string const& s) :
            std::runtime_error(s) {}
        UnitCellError(const char* s) :
            std::runtime_error(std::string(s)) {}
};

class InvalidInteger : public std::runtime_error {
   public:
      InvalidInteger(std::string const& s) :
         std::runtime_error(s) {}
};

class InvalidDecimal : public std::runtime_error {
   public:
      InvalidDecimal(std::string const& s) :
         std::runtime_error(s) {}
};

class StringBufferOverflow : public std::runtime_error {
   public:
      StringBufferOverflow(std::string const& s) :
         std::runtime_error(s) {}
};

class FileIOError : public std::runtime_error {
   public:
      FileIOError(std::string const& s) :
         std::runtime_error(s) {}
};

class InternalError : public std::runtime_error {
   public:
      InternalError(std::string const& s) :
         std::runtime_error(s) {}
};

};
#endif /* EXCEPTIONS_H */
