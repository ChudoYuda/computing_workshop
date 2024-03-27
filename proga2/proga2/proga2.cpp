
#include <string>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>


using namespace std;

ifstream f1("input.txt", ios::in);
ofstream f2("output.txt", ios::out);

class Equation
{
public:
    virtual double principalRoot() = 0;
    virtual double sum() = 0;
    virtual void print() = 0;
    virtual void print1() = 0;
   
    ~Equation() {
        }
};
class Linear : public Equation 
{
protected:
    double a;
    double b;

public:
    Linear ( double a1, double b1): a(a1), b(b1) { } 

    ~Linear() {
        }

    double principalRoot() 
    {
        double x;
        if (b == 0 || a == 0)
            x = 0;
        else
            x = (- b / a);
        return x;
    }
    
    double sum()
    {
        return this->principalRoot();
    }
    
    void print() 
    {
        double x = this->principalRoot();
        f2 << a << "x + " << b << " = 0\t" << x << "\n";
        std::cout << a << "x + " << b << " = 0\t" << x << "\n";
    }
    void print1()
    {
       
        f2 << "(" << a << "x + " << b << ")";
        std::cout << "(" << a << "x + " << b << ")";
    }
};

class Trigonometric : public Equation
{
protected:
    double k;
public:
    Trigonometric(double k1) : k(k1) { }

    ~Trigonometric() {
        }

    double principalRoot()
    {
        double x = (asin(k));
        return x;
    }

    double sum()
    {
        return this->principalRoot();
    }

    void print()
    {
        double x = this->principalRoot();
        f2 << "sin(x) - " << k << " = 0\t" << x << "\n";
        std::cout << "sin(x) - " << k << " = 0\t" << x << "\n";
    }
    void print1()
    {
        
        f2 << "(sin(x) - " << k << ")";
        std::cout << "(sin(x) - " << k << ")";
    }
};

class Disjunction: public Equation 
{
protected:
    Equation* equation1;
    Equation* equation2;

public:
    Disjunction (Equation* z1, Equation* z2) : equation1(z1), equation2(z2) { }

    double principalRoot()
    {
        double root1 = equation1->principalRoot();
        double root2 = equation2->principalRoot();

        if (abs(root1) < abs(root2))
             return root1;
        if (abs(root1) > abs(root2))
            return root2;
        else
            return abs(root1);
    }

    ~Disjunction() 
    {
        delete equation1;
        delete equation2;
    }

    double sum() 
    {
        return (equation1 -> sum() + equation2 -> sum());
    }

    void print()
    {
        
        equation1->print1();
        equation2->print1();
        double x = this->principalRoot();
        f2 << " = 0\t" << x << "\n";
        std::cout << " = 0\t" << x << "\n";
    }
    void print1()
    {
       
        equation1->print1();
        equation2->print1();
    }
};


int main()
{
    double N;
    double a;
    double b;
    double k;
    f1 >> N;
    char type;
    f1 >> type;
    Equation* ex1 = nullptr;
    Equation* ex2 = nullptr;
    switch (type)
    {
    case 'L':
        f1 >> a;
        f1 >> b;
        ex1 = new Linear(a, b);
        break;
    case 'T':
        f1 >> k;
        ex1 = new Trigonometric(k);
        break;
    default:
        f2 << "ERROR";
    }
    //ex1->print();
    for (int i = 0; i < N-1; i++)
    {
        f1 >> type;
        switch (type)
        {
        case 'L':
            f1 >> a;
            f1 >> b;
            ex2 = new Linear( a, b);
            break;
        case 'T':
            f1 >> k;
            ex2 = new Trigonometric(k);
            break;
        default:
            f2 << "ERROR";
        }
        //Equation* cop = ex1;
        //delete ex1;
        ex1 = new Disjunction(ex2, ex1); //не будет ли мусора в памяти?
    }
    f2 << fixed << setprecision(2) << (ex1->sum()) <<"\n";
    //ex1->print();
    delete ex1;
    
    return 0;
}

