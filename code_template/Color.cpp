#include <iomanip>
#include "Color.h"

Color::Color() {
    this->r = 0;
    this->g = 0;
    this->b = 0;
}

Color::Color(double r, double g, double b)
{
    this->r = r;
    this->g = g;
    this->b = b;
}

Color::Color(const Color &other)
{
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;
}

Color Color::operator+(const Color &other){
    Color result((this->r + other.r),(this->g + other.g),(this->b + other.b) );
    return result;
}

Color Color::operator-(const Color &other){
    Color result((this->r - other.r),(this->g - other.g),(this->b - other.b) );
    return result;
}

Color Color::operator*(const double &other){
    Color result((this->r * other),(this->g*other),(this->b*other));
    return result;
}

Color Color::operator/(const double &other){
    Color result((this->r / other),(this->g/other),(this->b/other));
    return result;
}

std::ostream &operator<<(std::ostream &os, const Color &c)
{
    os << std::fixed << std::setprecision(0) << "rgb(" << c.r << ", " << c.g << ", " << c.b << ")";
    return os;
}
