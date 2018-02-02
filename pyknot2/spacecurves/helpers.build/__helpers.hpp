#ifndef __NUITKA_CALLS_H__
#define __NUITKA_CALLS_H__

extern PyObject *CALL_FUNCTION_WITH_ARGS1( PyObject *called, PyObject *arg0 );
extern PyObject *CALL_FUNCTION_WITH_ARGS2( PyObject *called, PyObject *arg0, PyObject *arg1 );
extern PyObject *CALL_FUNCTION_WITH_ARGS3( PyObject *called, PyObject *arg0, PyObject *arg1, PyObject *arg2 );
extern PyObject *CALL_FUNCTION_WITH_ARGS4( PyObject *called, PyObject *arg0, PyObject *arg1, PyObject *arg2, PyObject *arg3 );
extern PyObject *CALL_FUNCTION_WITH_ARGS8( PyObject *called, PyObject *arg0, PyObject *arg1, PyObject *arg2, PyObject *arg3, PyObject *arg4, PyObject *arg5, PyObject *arg6, PyObject *arg7 );
#endif
