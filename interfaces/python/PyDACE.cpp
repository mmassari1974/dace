#include <dace/dace.h>
#include <dace/Monomial.h>
#include <dace/Interval.h>
#include <dace/compiledDA.h>
#include <dace/AlgebraicVector.h>
using namespace DACE;

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/stl_iterator.hpp>
using namespace boost::python;

/////////////////////////////////////////////////
//           CONVERTERS DEFINITION             //
/////////////////////////////////////////////////

// Function used to convert a C++ vector to a python list.
template <typename T>
boost::python::list vect2list(std::vector<T> InputVect)
{
  typename std::vector<T>::iterator iter;
  boost::python::list OutputList;

  for (iter = InputVect.begin(); iter != InputVect.end(); ++iter)
  {
    OutputList.append(*iter);
  }

  return OutputList;
};

template <typename T>
std::vector<T> list2vect(const boost::python::list &iterable)
{
  return std::vector<T>(boost::python::stl_input_iterator<T>(iterable),
                        boost::python::stl_input_iterator<T>());
}

// Functions used to convert Python lists to C++ vectors.
// More than one converter is needed because automatic
// conversion from a list to a specified C++ type is not
// possible. E.g. a list made of positive integers can
// be converted both in std::vector<int> and
// std::vector<unsigned int>. Same holds for floats and
// doubles.
// std::vector<int> list2vect_int(const boost::python::list &iterable)
// {
//   return std::vector<int>(boost::python::stl_input_iterator<int>(iterable),
//                           boost::python::stl_input_iterator<int>());
// }

// std::vector<unsigned int> list2vect_unsigned_int(const boost::python::list &iterable)
// {
//   return std::vector<unsigned int>(boost::python::stl_input_iterator<unsigned int>(iterable),
//                                    boost::python::stl_input_iterator<unsigned int>());
// }

// std::vector<float> list2vect_float(const boost::python::list &iterable)
// {
//   return std::vector<float>(boost::python::stl_input_iterator<float>(iterable),
//                             boost::python::stl_input_iterator<float>());
// }

// std::vector<double> list2vect_double(const boost::python::list &iterable)
// {
//   return std::vector<double>(boost::python::stl_input_iterator<double>(iterable),
//                              boost::python::stl_input_iterator<double>());
// }

// std::vector<std::string> list2vect_string(const boost::python::list &iterable)
// {
//   return std::vector<std::string>(boost::python::stl_input_iterator<std::string>(iterable),
//                                   boost::python::stl_input_iterator<std::string>());
// }

// std::vector<DA> list2vect_DA(const boost::python::list &iterable)
// {
//   return std::vector<DA>(boost::python::stl_input_iterator<DA>(iterable),
//                          boost::python::stl_input_iterator<DA>());
// }

// Function used to print C++ vectors in Python
// template <typename T>
// void printVect(std::vector<T> InputVect)
// {
//   std::cout << "[";

//   for (int i = 0; i != InputVect.size() - 1; ++i)
//   {
//     std::cout << InputVect[i] << ", ";
//   }

//   std::cout << InputVect[InputVect.size() - 1] << "] \n";
// };

/////////////////////////////////////////////////
//         WRAPPING WITH BOOST.PYTHON          //
/////////////////////////////////////////////////

class PyDA : public DA
{
  private:
  protected:
  public:

  boost::python::list linear() const
  {
    return vect2list<double>(DA::linear());
  }

  boost::python::list gradient() const
  {
    return vect2list<PyDA>(DA::gradient());
  }

  double getCoefficient(const boost::python::list &jj) const
  {
    return DA::getCoefficient(list2vect<unsigned int>(jj));
  }                //!< Get specific coefficient
  void setCoefficient(const boost::python::list &jj, const double coeff)
  {
    DA::setCoefficient(list2vect<unsigned int>(jj) , coeff);
  }

};


BOOST_PYTHON_MODULE(PyDACE)
{
  /////////////////////////////////////////////////
  //  WRAPPING OF THE AlgebraicVector, Monomial  //
  //  Interval AND compiledDA CLASSES            //
  /////////////////////////////////////////////////

  // Since AlgebraicVector is a templated class, for each type
  // a wrap is carried out
  class_<AlgebraicVector<int>>("AlgebraicVector_int", init<>());
  class_<AlgebraicVector<unsigned int>>("AlgebraicVector_unsigned_int", init<>());
  class_<AlgebraicVector<float>>("AlgebraicVector_float", init<>());
  class_<AlgebraicVector<double>>("AlgebraicVector_double", init<>());
  class_<AlgebraicVector<DA>>("AlgebraicVector_DA", init<>());
  class_<AlgebraicVector<Monomial>>("AlgebraicVector_Monomial", init<>());

  class_<Monomial>("Monomial", init<>())
      .def("order", &Monomial::order)
      .def("toString", &Monomial::toString);

  class_<Interval>("Interval", init<>())
      .def_readonly("m_lb", &Interval::m_lb)
      .def_readonly("m_ub", &Interval::m_ub);

  class_<compiledDA>("compiledDA", init<const compiledDA>())
      .def(init<const DA>())
      .def(init<const std::vector<DA>>());

  /////////////////////////////////////////////////
  //         WRAPPING OF THE CONVERTERS          //
  /////////////////////////////////////////////////

  // All overloaded methods are wrapped
  // def("vect2list", vect2list<int>);
  // def("vect2list", vect2list<unsigned int>);
  // def("vect2list", vect2list<float>);
  // def("vect2list", vect2list<double>);
  // def("vect2list", vect2list<std::string>);
  // def("vect2list", vect2list<DA>);

  // def("list2vect_int", list2vect_int);
  // def("list2vect_unsigned_int", list2vect_unsigned_int);
  // def("list2vect_float", list2vect_float);
  // def("list2vect_double", list2vect_double);
  // def("list2vect_string", list2vect_string);
  // def("list2vect_DA", list2vect_DA);

  // def("printVect", printVect<int>);
  // def("printVect", printVect<float>);
  // def("printVect", printVect<unsigned int>);
  // def("printVect", printVect<double>);
  // def("printVect", printVect<std::string>);

  /////////////////////////////////////////////////
  //    WRAPPING OF THE std::vector<> CLASSES    //
  /////////////////////////////////////////////////

  // As for the AlgebraicVector, for each type of
  // std::vector a wrap is carried out
  // class_<std::vector<int>>("stdVector_int")
  //     .def(vector_indexing_suite<std::vector<int>>());

  // class_<std::vector<unsigned int>>("stdVector_unsigned_int")
  //     .def(vector_indexing_suite<std::vector<unsigned int>>());

  // class_<std::vector<float>>("stdVector_float")
  //     .def(vector_indexing_suite<std::vector<float>>());

  // class_<std::vector<double>>("stdVector_double")
  //     .def(vector_indexing_suite<std::vector<double>>());

  // class_<std::vector<std::string>>("stdVector_string")
  //     .def(vector_indexing_suite<std::vector<std::string>>());

  // Vectors of Monomial and DA are not accepted,
  // because of the reason explained here:
  // http://stackoverflow.com/questions/10680691/why-do-i-need-comparison-operators-in-boost-python-vector-indexing-suite
  // To create these a little modification in the
  // Monomial.h and DA.h codes are needed.
  //
  // class_< std::vector<Monomial> >("stdVector_Monomial")
  // .def(vector_indexing_suite<std::vector<Monomial> >());
  //
  // class_< std::vector<DA> >("stdVector_DA")
  // .def(vector_indexing_suite<std::vector<DA> >());

  /////////////////////////////////////////////////
  // DEFINITIONS OF FUNCTION METHODS OVERLOADING //
  /////////////////////////////////////////////////

  // Overloading made to consider also the
  // user-friendly implementation of the functions.
  // "nm" and "m" stand for non/member function.

  // double (DA::*m_cons)() const = &DA::cons;
  // double (*nm_cons)(const DA &) = cons;

  // // Missing "linear"
  // // Missing "gradient"

  // DA (DA::*m_derivSingle)
  // (const unsigned int) const = &DA::deriv;
  // DA (*nm_derivSingle)
  // (const DA &, const unsigned int) = deriv;

  // DA (DA::*m_derivMultiple)
  // (const std::vector<unsigned int>) const = &DA::deriv;
  // DA (*nm_derivMultiple)
  // (const DA &, const std::vector<unsigned int>) = deriv;

  // DA (DA::*m_integSingle)
  // (const unsigned int) const = &DA::integ;
  // DA (*nm_integSingle)
  // (const DA &, const unsigned int) = integ;

  // DA (DA::*m_integMultiple)
  // (const std::vector<unsigned int>) const = &DA::integ;
  // DA (*nm_integMultiple)
  // (const DA &, const std::vector<unsigned int>) = integ;

  // DA (DA::*m_trim)
  // (const unsigned int, const unsigned int) const = &DA::trim;
  // DA (*nm_trim)
  // (const DA &, const unsigned int, const unsigned int) = trim;

  // DA (DA::*m_trunc)
  // () const = &DA::trunc;
  // DA (*nm_trunc)
  // (const DA &) = trunc;

  // DA (DA::*m_round)
  // () const = &DA::round;
  // DA (*nm_round)
  // (const DA &) = round;

  // // ERROR!
  // // Inside the DA class there exists the
  // // member function "mod", while the nonmember
  // // function is called "fmod". Looks like fmod
  // // has its own declaration in DA.h but not a
  // // definition. Probably it should be named "mod"
  // // instead of "fmod".
  // //
  // DA (DA::*m_mod)
  // (const double) const = &DA::mod;
  // // DA (*nm_mod)(const DA&, const double) = fmod;

  // DA (DA::*m_pow)
  // (const int) const = &DA::pow;
  // DA (*nm_pow)
  // (const DA &, const int) = pow;

  // DA (DA::*m_root)
  // (const int) const = &DA::root;
  // DA (*nm_root)
  // (const DA &, const int) = root;

  // DA (DA::*m_minv)
  // () const = &DA::minv;
  // DA (*nm_minv)
  // (const DA &) = minv;

  // DA (DA::*m_sqr)
  // () const = &DA::sqr;
  // DA (*nm_sqr)
  // (const DA &) = sqr;

  // DA (DA::*m_sqrt)
  // () const = &DA::sqrt;
  // DA (*nm_sqrt)
  // (const DA &) = sqrt;

  // DA (DA::*m_isrt)
  // () const = &DA::isrt;
  // DA (*nm_isrt)
  // (const DA &) = isrt;

  // DA (DA::*m_exp)
  // () const = &DA::exp;
  // DA (*nm_exp)
  // (const DA &) = exp;

  // DA (DA::*m_log)
  // () const = &DA::log;
  // DA (*nm_log)
  // (const DA &) = log;

  // DA (DA::*m_logb)
  // (const double) const = &DA::logb;
  // DA (*nm_logb)
  // (const DA &, const double) = logb;

  // DA (DA::*m_sin)
  // () const = &DA::sin;
  // DA (*nm_sin)
  // (const DA &) = sin;

  // DA (DA::*m_cos)
  // () const = &DA::cos;
  // DA (*nm_cos)
  // (const DA &) = cos;

  // DA (DA::*m_tan)
  // () const = &DA::tan;
  // DA (*nm_tan)
  // (const DA &) = tan;

  // DA (DA::*m_asin)
  // () const = &DA::asin;
  // DA (*nm_asin)
  // (const DA &) = asin;

  // DA (DA::*m_acos)
  // () const = &DA::acos;
  // DA (*nm_acos)
  // (const DA &) = acos;

  // DA (DA::*m_atan)
  // () const = &DA::atan;
  // DA (*nm_atan)
  // (const DA &) = atan;

  // DA (DA::*m_atan2)
  // (const DA &) const = &DA::atan2;
  // DA (*nm_atan2)
  // (const DA &, const DA &) = atan2;

  // DA (DA::*m_sinh)
  // () const = &DA::sinh;
  // DA (*nm_sinh)
  // (const DA &) = sinh;

  // DA (DA::*m_cosh)
  // () const = &DA::cosh;
  // DA (*nm_cosh)
  // (const DA &) = cosh;

  // DA (DA::*m_tanh)
  // () const = &DA::tanh;
  // DA (*nm_tanh)
  // (const DA &) = tanh;

  // DA (DA::*m_asinh)
  // () const = &DA::asinh;
  // DA (*nm_asinh)
  // (const DA &) = asinh;

  // DA (DA::*m_acosh)
  // () const = &DA::acosh;
  // DA (*nm_acosh)
  // (const DA &) = acosh;

  // DA (DA::*m_atanh)
  // () const = &DA::atanh;
  // DA (*nm_atanh)
  // (const DA &) = atanh;

  // unsigned int (DA::*m_size)() const = &DA::size;
  // unsigned int (*nm_size)(const DA &) = size;

  // double (DA::*m_abs)() const = &DA::abs;
  // double (*nm_abs)(const DA &) = abs;

  // double (DA::*m_norm)(unsigned int) const = &DA::norm;
  // double (*nm_norm)(const DA &, unsigned int) = norm;

  // std::vector<double> (DA::*m_orderNorm)(unsigned int, unsigned int) const = &DA::orderNorm;
  // std::vector<double> (*nm_orderNorm)(const DA &, unsigned int, unsigned int) = orderNorm;

  // std::vector<double> (DA::*m_estimNorm)(unsigned int, unsigned int, unsigned int) const = &DA::estimNorm;
  // std::vector<double> (*nm_estimNorm)(const DA &, unsigned int, unsigned int, unsigned int) = estimNorm;

  // Interval (DA::*m_bound)() const = &DA::bound;
  // Interval (*nm_bound)(const DA &) = bound;

  // double (DA::*m_convRadius)(const double, const unsigned int) const = &DA::convRadius;
  // double (*nm_convRadius)(const DA &, const double, const unsigned int) = convRadius;

  // // Definition of overloaded different methods
  // // for the Template functions. Floats and doubles
  // // are considered by Python as equal; nonetheless
  // // two different methods are provided hereinafter.
  // // (May be avoided?)
  // // Moreover, is it useful to define arrays for Python?
  // // (Also the evalArray functions may be avoided)

  // int (DA::*m_evalVector_int)(const std::vector<int> &) const = &DA::eval;
  // int (*nm_evalVector_int)(const DA &, const std::vector<int> &) = eval;

  // unsigned int (DA::*m_evalVector_unsigned_int)(const std::vector<unsigned int> &) const = &DA::eval;
  // unsigned int (*nm_evalVector_unsigned_int)(const DA &, const std::vector<unsigned int> &) = eval;

  // float (DA::*m_evalVector_float)(const std::vector<float> &) const = &DA::eval;
  // float (*nm_evalVector_float)(const DA &, const std::vector<float> &) = eval;

  // double (DA::*m_evalVector_double)(const std::vector<double> &) const = &DA::eval;
  // double (*nm_evalVector_double)(const DA &, const std::vector<double> &) = eval;

  // int (DA::*m_evalArray_int)(const int[], const unsigned int) const = &DA::eval;
  // int (*nm_evalArray_int)(const DA &, const int[], const unsigned int) = eval;

  // unsigned int (DA::*m_evalArray_unsigned_int)(const unsigned int[], const unsigned int) const = &DA::eval;
  // unsigned int (*nm_evalArray_unsigned_int)(const DA &, const unsigned int[], const unsigned int) = eval;

  // float (DA::*m_evalArray_float)(const float[], const unsigned int) const = &DA::eval;
  // float (*nm_evalArray_float)(const DA &, const float[], const unsigned int) = eval;

  // double (DA::*m_evalArray_double)(const double[], const unsigned int) const = &DA::eval;
  // double (*nm_evalArray_double)(const DA &, const double[], const unsigned int) = eval;

  // int (DA::*m_evalScalar_int)(const int &) const = &DA::evalScalar;
  // int (*nm_evalScalar_int)(const DA &, const int &) = evalScalar;

  // unsigned int (DA::*m_evalScalar_unsigned_int)(const unsigned int &) const = &DA::evalScalar;
  // unsigned int (*nm_evalScalar_unsigned_int)(const DA &, const unsigned int &) = evalScalar;

  // float (DA::*m_evalScalar_float)(const float &) const = &DA::evalScalar;
  // float (*nm_evalScalar_float)(const DA &, const float &) = evalScalar;

  // double (DA::*m_evalScalar_double)(const double &) const = &DA::evalScalar;
  // double (*nm_evalScalar_double)(const DA &, const double &) = evalScalar;

  // compiledDA (DA::*m_compile)() const = &DA::compile;
  // compiledDA (*nm_compile)(const DA &) = compile;

  // DA (DA::*m_plug)
  // (unsigned int, double) const = &DA::plug;
  // DA (*nm_plug)
  // (const DA &, unsigned int, double) = plug;

  // std::string (DA::*m_toString)() const = &DA::toString;
  // std::string (*nm_toString)(const DA &) = toString;

  // Missing fromString, both for the single string
  // and the std::vect<string> inputs.
  //
  // DA (DA::*fromStringSingle)(const std::string&) = &DA::fromString;
  // DA (DA::*fromStringMultiple)(const std::vector<std::string>&) = &DA::fromString;

  /////////////////////////////////////////////////
  //          WRAPPING OF THE DA CLASS           //
  /////////////////////////////////////////////////

  /********************************************************************************
    *     Constructors & Destructors
    *********************************************************************************/
  // Attention!
  // The order in which the constructors are defined
  // hereinafter matters. They're considered starting
  // from the last one moving to the first. Hence,
  // keep this order and don't change it!
  class_<PyDA>("DA", init<>())
      .def(init<const double>())
      .def(init<const PyDA>())
      .def(init<const unsigned int, optional<const double>>())
      .def(init<const int, optional<const double>>())

      /********************************************************************************
    *     DACE Setup
    *********************************************************************************/
      .def("init", &DA::init)
      .def("isInitialized", &DA::isInitialized)

      // Error in testing/don't know how to test
      // the "version" function.
      .def("version", &DA::version)

      .def("checkVersion", &DA::checkVersion)
      .def("setEps", &DA::setEps)
      .def("getEps", &DA::getEps)
      .def("getEpsMac", &DA::getEpsMac)
      .def("getMaxOrder", &DA::getMaxOrder)
      .def("getMaxVariables", &DA::getMaxVariables)
      .def("getMaxMonomials", &DA::getMaxMonomials)
      .def("getTO", &DA::getTO)
      .def("setTO", &DA::setTO)
      .def("pushTO", &DA::pushTO)
      .def("popTO", &DA::popTO)

      .staticmethod("init")
      .staticmethod("isInitialized")
      .staticmethod("version")
      .staticmethod("checkVersion")
      .staticmethod("setEps")
      .staticmethod("getEps")
      .staticmethod("getEpsMac")
      .staticmethod("getMaxOrder")
      .staticmethod("getMaxVariables")
      .staticmethod("getMaxMonomials")
      .staticmethod("getTO")
      .staticmethod("setTO")
      .staticmethod("pushTO")
      .staticmethod("popTO")

      /********************************************************************************
    *     Coefficient access and extraction routines
    *********************************************************************************/
      .def("cons", m_cons)

      // Watch out. AlgebraicVector here
      // used for "linear" and "gradient".
      .def("linear", &DA::linear)
      .def("gradient", &DA::gradient)

      .def("getCoefficient", &DA::getCoefficient)
      .def("setCoefficient", &DA::setCoefficient)

      .def("getMonomial", static_cast<Monomial (DA::*)(const unsigned int) const>(&DA::getMonomial))
      .def("getMonomials", &DA::getMonomials)

      /********************************************************************************
    *     Assignments
    *********************************************************************************/

      .def(self += self)
      .def(self += double())

      .def(self -= self)
      .def(self -= double())

      .def(self *= self)
      .def(self *= double())

      .def(self /= self)
      .def(self /= double())

      /********************************************************************************
    *     Algebraic operations
    *********************************************************************************/

      .def(-self)

      .def(self + self)
      .def(self + double())
      .def(double() + self)

      .def(self - self)
      .def(self - double())
      .def(double() - self)

      .def(self * self)
      .def(self * double())
      .def(double() * self)

      .def(self / self)
      .def(self / double())
      .def(double() / self)

      /********************************************************************************
    *     Math routines
    *********************************************************************************/

      .def("deriv", m_derivSingle)
      .def("deriv", m_derivMultiple)
      .def("integ", m_integSingle)
      .def("integ", m_integMultiple)

      .def("trim", m_trim, (arg("max") = DA::getMaxOrder()))
      .def("trunc", m_trunc)
      .def("round", m_round)
      .def("mod", m_mod)
      .def("pow", m_pow)
      .def("root", m_root, (arg("p") = 2))
      .def("minv", m_minv)
      .def("sqr", m_sqr)
      .def("sqrt", m_sqrt)
      .def("isrt", m_isrt)
      .def("exp", m_exp)
      .def("log", m_log)
      .def("logb", m_logb, (arg("b") = 10.0))
      .def("sin", m_sin)
      .def("cos", m_cos)
      .def("tan", m_tan)
      .def("asin", m_asin)
      .def("acos", m_acos)
      .def("atan", m_atan)
      .def("atan2", m_atan2)
      .def("sinh", m_sinh)
      .def("cosh", m_cosh)
      .def("tanh", m_tanh)
      .def("asinh", m_asinh)
      .def("acosh", m_acosh)
      .def("atanh", m_atanh)

      /********************************************************************************
    *    Norm and estimation routines
    *********************************************************************************/

      .def("size", m_size)
      .def("abs", m_abs)
      .def("norm", m_norm, (arg("type") = 0))
      .def("orderNorm", m_orderNorm, (arg("var") = 0, arg("type") = 0))
      .def("estimNorm", m_estimNorm, (arg("var") = 0, arg("type") = 0, arg("nc") = DA::getMaxOrder()))

      .def("bound", m_bound)

      .def("convRadius", m_convRadius, (arg("type") = 1))

      /********************************************************************************
    *     DACE polynomial evaluation routines
    *********************************************************************************/

      // Watch out. "eval" and "evalScalar" use
      // templates as input and return types.
      .def("eval", m_evalVector_int)
      .def("eval", m_evalVector_unsigned_int)
      .def("eval", m_evalVector_float)
      .def("eval", m_evalVector_double)

      .def("eval", m_evalArray_int)
      .def("eval", m_evalArray_unsigned_int)
      .def("eval", m_evalArray_float)
      .def("eval", m_evalArray_double)

      .def("evalScalar", m_evalScalar_int)
      .def("evalScalar", m_evalScalar_unsigned_int)
      .def("evalScalar", m_evalScalar_float)
      .def("evalScalar", m_evalScalar_double)

      .def("compile", m_compile)

      .def("plug", m_plug, (arg("val") = 0.0))

      /********************************************************************************
    *     DACE input/output routines
    *********************************************************************************/

      .def("toString", m_toString)

      // No need to wrap also the istream and
      // ostream operators.

      /********************************************************************************
    *     Static factory routines
    *********************************************************************************/

      .def("random", &DA::random)
      .def("identity", &DA::identity)

      // Error due to the C++ compiler and the
      // conversion from std::string.
      // After having solved this problem, take
      // of the overloading, since both a 1D and
      // a multidimensional version of this
      // function exist.
      .def("fromString", static_cast<DA (*)(const std::string &)>(&DA::fromString))

      .staticmethod("random")
      .staticmethod("identity")

      /********************************************************************************
    *     DACE various routines
    *********************************************************************************/

      // How can I test this function?
      .def("memdump", &DA::memdump)
      .staticmethod("memdump");

  /********************************************************************************
    *     DACE non-member functions
    *********************************************************************************/

  def("cons", nm_cons);
  def("deriv", nm_derivSingle);
  def("deriv", nm_derivMultiple);
  def("integ", nm_integSingle);
  def("integ", nm_integMultiple);
  def("trim", nm_trim, (arg("max") = DA::getMaxOrder()));
  def("trunc", nm_trunc);
  def("round", nm_round);

  // Error mentioned above
  // def("mod", nm_mod);

  def("pow", nm_pow);
  def("root", nm_root, (arg("p") = 2));
  def("minv", nm_minv);
  def("sqr", nm_sqr);
  def("sqrt", nm_sqrt);
  def("isrt", nm_isrt);
  def("exp", nm_exp);
  def("log", nm_log);
  def("logb", nm_logb, (arg("b") = 10.0));
  def("sin", nm_sin);
  def("cos", nm_cos);
  def("tan", nm_tan);
  def("asin", nm_asin);
  def("acos", nm_acos);
  def("atan", nm_atan);
  def("atan2", nm_atan2);
  def("sinh", nm_sinh);
  def("cosh", nm_cosh);
  def("tanh", nm_tanh);
  def("asinh", nm_asinh);
  def("acosh", nm_acosh);
  def("atanh", nm_atanh);

  def("size", nm_size);
  def("abs", nm_abs);
  def("norm", nm_norm, (arg("type") = 0));
  def("orderNorm", nm_orderNorm, (arg("var") = 0, arg("type") = 0));
  def("estimNorm", nm_estimNorm, (arg("var") = 0, arg("type") = 0, arg("nc") = DA::getMaxOrder()));
  def("bound", nm_bound);
  def("convRadius", nm_convRadius, (arg("type") = 1));

  def("eval", nm_evalVector_int);
  def("eval", nm_evalVector_unsigned_int);
  def("eval", nm_evalVector_float);
  def("eval", nm_evalVector_double);

  def("eval", nm_evalArray_int);
  def("eval", nm_evalArray_unsigned_int);
  def("eval", nm_evalArray_float);
  def("eval", nm_evalArray_double);

  def("evalScalar", nm_evalScalar_int);
  def("evalScalar", nm_evalScalar_unsigned_int);
  def("evalScalar", nm_evalScalar_float);
  def("evalScalar", nm_evalScalar_double);

  def("compile", nm_compile);
  def("plug", nm_plug, (arg("val") = 0.0));
  def("toString", nm_toString);
}