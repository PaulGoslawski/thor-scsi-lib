#ifndef _TPS_SS_VECT_H_
#define _TPS_SS_VECT_H_

#include <tps/config.h>
#include <tps/enums.h>
#include <tps/forward_decl.h>
#include <tps/tps_type.h>
#include <string>
#include <vector>


// pre-declare template friend functions: f<>()
template<typename T>
ss_vect<T> operator+(const ss_vect<T> &);
template<typename T>
ss_vect<T> operator-(const ss_vect<T> &);


template<typename T> class ss_vect {
 public:
  typedef T value_type;

  ss_vect(void);
  ss_vect(const double x, const double px, const double y, const double py,
	  const double ct, const double delta);
// Let's the compiler synthetize the copy constructor
//  ss_vect(const T &a) { }
//  ss_vect(const ss_vect<T> &a) { }
  template<typename U>
    ss_vect(const U &);
  template<typename U>
    ss_vect(const ss_vect<U> &);

  void print(const std::string &);

  ss_vect<double> cst(void) const;
  T& operator[](const int i) { return ss[i]; }
  const T& operator[](const int i) const { return ss[i]; }

  ss_vect<T>& operator*=(const double);
  ss_vect<T>& operator*=(const tps &);

  ss_vect<T>& operator=(const ss_vect<T> &);
  ss_vect<T>& operator+=(const ss_vect<T> &);
  ss_vect<T>& operator-=(const ss_vect<T> &);

  friend ss_vect<T> operator+<>(const ss_vect<T> &);
  friend ss_vect<T> operator-<>(const ss_vect<T> &);

  friend ss_vect<double> operator+(const ss_vect<double> &, const ss_vect<double> &);
  friend ss_vect<tps> operator+(const ss_vect<tps> &, const ss_vect<double> &);
  friend ss_vect<tps> operator+(const ss_vect<double> &, const ss_vect<tps> &);
  friend ss_vect<tps> operator+(const ss_vect<tps> &, const ss_vect<tps> &);

  friend ss_vect<double> operator-(const ss_vect<double> &, const ss_vect<double> &);
  friend ss_vect<tps> operator-(const ss_vect<tps> &, const ss_vect<double> &);
  friend ss_vect<tps> operator-(const ss_vect<double> &, const ss_vect<tps> &);
  friend ss_vect<tps> operator-(const ss_vect<tps> &, const ss_vect<tps> &);

//  friend ss_vect<double> operator*(const ss_vect<tps> &,
//				   const ss_vect<double> &);
  // R(nd2, nv) = P(nd2, nd2)*Q(nd2, nv)
  friend ss_vect<double> operator*(const double, const ss_vect<double> &);
  friend ss_vect<double> operator*(const ss_vect<double> &, const double);
  friend ss_vect<tps> operator*(const ss_vect<tps> &, const ss_vect<tps> &);
  friend ss_vect<tps> operator*(const double, const ss_vect<tps> &);
  friend ss_vect<tps> operator*(const ss_vect<tps> &, const double);
  // R(nv, nv) = P(nv, nv)*Q(nv, nv)
  friend ss_vect<tps> CCT(const ss_vect<tps> &, const ss_vect<tps> &);
  friend ss_vect<tps> MTREE(const ss_vect<tps> &);
  friend ss_vect<double> PPUSH(const ss_vect<tps> &, ss_vect<double> &);
  friend tps operator*(const tps &, const ss_vect<tps> &);

  template<typename CharT, class Traits>
    friend std::basic_istream<CharT, Traits>&
    operator>>(std::basic_istream<CharT, Traits> &, ss_vect<tps> &);

  template<typename CharT, class Traits>
    friend std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits> &, const ss_vect<T> &);

  ss_vect<T> zero(void);
  ss_vect<T> identity(void);

  friend ss_vect<tps> FExpo(const tps &, const ss_vect<tps> &,
			    const int, const int, const int);
  friend ss_vect<tps> LieExp(const tps &, const ss_vect<tps> &);
  friend ss_vect<tps> LieFlo(const ss_vect<tps> &, const ss_vect<tps> &);
  // Q(nv, nv) = P(nd2, nd2)^-1
  friend ss_vect<tps> Inv(const ss_vect<tps> &);
  // Q(nv, nv) = P(nv, nv)^-1
  friend ss_vect<tps> Inv_Ext(const ss_vect<tps> &);
  friend ss_vect<tps> PInv(const ss_vect<tps> &, const long int []);
  friend void GoFix(const ss_vect<tps> &, ss_vect<tps> &,
		    ss_vect<tps> &, const int);
  friend MNF_struct MapNorm(const ss_vect<tps> &, const int);
  friend ss_vect<tps> MapNormF(const ss_vect<tps> &, ss_vect<tps> &,
			       ss_vect<tps> &, ss_vect<tps> &,
			       ss_vect<tps> &, const int, const int);
  friend void dHdJ(const tps &, ss_vect<tps> &);
  friend void CtoR(const tps &, tps &, tps &);
  friend tps RtoC(const tps &, const tps &);
  friend tps LieFact_DF(const ss_vect<tps> &, ss_vect<tps> &);
  friend tps LieFact(const ss_vect<tps> &);
  friend ss_vect<tps> FlowFact(const ss_vect<tps> &);
  friend tps Intd(const ss_vect<tps> &, const double);
  friend ss_vect<tps> Difd(const tps &, const double);
  friend ss_vect<tps> Taked(const ss_vect<tps> &, const int);


 private:
  // (Note, e.g. spin components should be added here)
  std::vector<T> ss;
};


// Class for single particle phase space dynamics



template<>
inline ss_vect<double> ss_vect<tps>::cst(void) const
{
  int             i;
  ss_vect<double> x;

  for (i = 0; i < ps_dim; i++)
    x[i] = (*this)[i].cst();
  return x;
}


// partial specialization
template<>
class is_double< ss_vect<double> > {
 public:
  static inline ss_vect<double> cst(const ss_vect<double> &x) { return x; }
  static inline ss_vect<double> ps(const ss_vect<tps> &x) { return x.cst(); }
};

// partial specialization
template<>
class is_double< ss_vect<tps> > {
 public:
  static inline ss_vect<double> cst(const ss_vect<tps> &x) { return x.cst(); }
  static inline ss_vect<tps> ps(const ss_vect<tps> &x) { return x; }
};
template<typename T>
inline ss_vect<T>::ss_vect(void) { ss.resize(ps_dim); }

template<typename T>
inline ss_vect<T>::ss_vect(const double x, const double px, const double y,
			   const double py, const double delta, const double ct)
{ ss = {x, px, y, py, delta, ct}; }

template<typename T>
template<typename U>
inline ss_vect<T>::ss_vect(const U &a) : ss(a) { }

template<typename T>
template<typename U>
inline ss_vect<T>::ss_vect(const ss_vect<U> &a)
{ for (int i = 0; i < ps_dim; i++) ss.push_back(a[i]); }



double xabs(long n, ss_vect<double> &x);

#endif /* _TPS_SS_VECT_H_ */