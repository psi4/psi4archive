#ifndef _psi_include_exception_h
#define _psi_include_exception_h

#include <exception>

namespace psi {

#define CHARARR_SIZE 100

#define PSIEXCEPTION(message) PsiException(message, __FILE__, __LINE__)

class PsiException : public std::runtime_error {
    private:
        std::string msg_;
        std::string file_;
        int line_;

    public:
        PsiException(
            std::string message,
            std::string file,
            int line
        ) throw () : std::runtime_error(message)
        {
            msg_ = message;
            file_ = file;
            line_ = line;
        }

        ~PsiException() throw ()
        {
        }

        const char* what() const throw ()
        {
            return msg_.c_str();
        }

        const char* file() const throw()
        {
            return file_._c_str();
        }

        int line() const throw()
        {
            return line_;
        }


};

class SanityCheckError : public PsiException {

    public:
        SanityCheckError(
            const char* message,
            const char* file,
            int line
            ) : PsiException(message, file, line)
        {
        }
};

template <
    class T
>
class LimitExceeded : public PsiException {

    private:
        T maxval_;
        T errorval_;

    public:
        LimitExceeded(
            const char* msg,
            T maxval,
            T errorval,
            const char* file,
            int line) : 
           PsiException(msg, file, line),
           maxval_(maxval), errorval_(errorval)
        {
        }

        T max_value(){return maxval_;}

        T actual_value(){return errorval_;}
};

class StepSizeError : public LimitExceeded<double> {

    typedef LimitExceeded<double> ParentClass;

    public:
        StepSizeError(
            const char* msg,
            double max,
            double actual,
            const char* file,
            int line)
            : ParentClass(msg, max, actual, file, line)
        {
        }

};


class MaxIterationsExceeded : public LimitExceeded<int> {

    typedef LimitExceeded<int> ParentClass;

    public:
        MaxIterationsExceeded(
            const char* msg,
            int max,
            const char* file,
            int line) :
            ParentClass(msg, max, max + 1, file, line)
        {
        }
};

class ConvergenceError : public MaxIterationsExceeded {

    public:
        ConvergenceError(
            const char* msg,
            int max,
            double desired_accuracy,
            double actual_accuracy,
            const char* file,
            int line) : 
            MaxIterationsExceeded(msg, max, file, line), desired_acc_(desired_accuracy), actual_acc_(actual_accuracy)
        {
        }

        double desired_accuracy(){return desired_acc_;}

        double actual_accuracy(){return actual_acc_;}

    private:
        double desired_acc_;
        double actual_acc_;


};

} //end namespace psi exception

#endif

