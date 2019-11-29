#[macro_use]
extern crate cpython;

use cpython::{Python, PyResult};


py_module_initializer!(libmyrustlib, initlibmyrustlib, PyInit_myrustlib, |py, m | {
    try!(m.add(py, "__doc__", "This module is implemented in Rust"));
    try!(m.add(py, "count_doubles", py_fn!(py, count_doubles(val: &str))));
    Ok(())
});