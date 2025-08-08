#include <iostream>
#include "cpplib.h"
#include "cxxbridge-cpp/lib.h"
#include "rust/cxx.h"

RsImage read_image(rust::Str path) {
    std::cout << "read_image called" << std::endl;
    std::cout << path << std::endl;
    Rgba c = { 1.0, 2.0, 3.0, 4.0};
    RsImage v = { 1, 1, c};
    return v;
}
void write_image(::rust::Str path, ::RsImage const & image) {

}
