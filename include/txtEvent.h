#ifndef TXTEVENT_H
#define TXTEVENT_H
#include <string>
#include <list>
#include <pputil.h>
#include <Vec.h>
#include <lorentz.h>
#include <particle.h>
    extern particleDataTable PDGtable;
    class particle;
    class txtEvent {
        protected:
      int _run;
      int _event;
      int _code;
      double _weight;
            std::list<particle> _initial;
            std::list<particle> _final;
            particle* _beam;
            particle* _target;
        public:
    txtEvent();
    txtEvent(const txtEvent&);
    ~txtEvent();
    txtEvent& operator=(const txtEvent& e);
    txtEvent& addfinal(const particle&);
    txtEvent& addinitial(const particle& p);
    txtEvent& erase();
    txtEvent& beam(const particle&);
    txtEvent& target(const particle&);
    int OK(double epsilon) const;
    particle beam() const;
    particle target() const;
    fourVec getPartPFinal(std::string name,int charge,int index,int debug=0) const;
    fourVec getPartPInitial(std::string name,int charge,int index) const;
    int f_charge() const;
    std::list<particle> f_mesons() const;
    std::list<particle> f_baryons() const;
    std::list<particle> f_particles() const;
    particle f_particle(const std::string& name, int charge, int index) const;
    int i_charge() const;
    std::list<particle> i_mesons() const;
    std::list<particle> i_baryons() const;
    std::list<particle> i_particles() const;
    particle i_particle(const std::string& name, int charge, int index) const;
    threeVec mesonPlane() const;
    threeVec baryonPlane() const;
    void print() const;
    friend std::istream& operator>>(std::istream& is, txtEvent& e);
    friend std::ostream& operator<<(std::ostream& os, txtEvent& e);
    friend txtEvent operator*(const lorentzTransform&,const txtEvent&);
    int run() const {return(_run);}
    int event() const {return(_event);}
    int run(int i) {return(_run = i);}
    int event(int i) {return(_event = i);}
    double weight() const {return(_weight);}
    int code() const {return(_code);}
    double weight(double w) {return(_weight = w);}
    int code(int i) {return(_code = i);}
    };
#endif
