//
// Copyright (C) 2011 Mateusz Loskot <mateusz@loskot.net>
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// Blog article: http://mateusz.loskot.net/?p=2819
#pragma once
#include <functional>
#include <iostream>
#include <string>
#include <Python.h>

namespace emb
{

typedef std::function<void(std::string)> stderr_write_type;

struct Stderr
{
  PyObject_HEAD
  stderr_write_type write;
};

typedef std::function<void(std::string)> stdout_write_type;

struct Stdout
{
  PyObject_HEAD
  stdout_write_type write;
};

PyObject* Stderr_write(PyObject* self, PyObject* args);
PyObject* Stderr_flush(PyObject* self, PyObject* args);
extern PyMethodDef Stderr_methods[];

PyObject* Stdout_write(PyObject* self, PyObject* args);
PyObject* Stdout_flush(PyObject* self, PyObject* args);
extern PyMethodDef Stdout_methods[];

extern PyModuleDef embmodule;

// Internal state
extern PyObject* g_stderr;
extern PyObject* g_stderr_saved;
extern PyObject* g_stdout;
extern PyObject* g_stdout_saved;

PyMODINIT_FUNC PyInit_emb(void);

void set_stderr(stderr_write_type write);
void reset_stderr();

void set_stdout(stdout_write_type write);
void reset_stdout();

} // namespace emb

