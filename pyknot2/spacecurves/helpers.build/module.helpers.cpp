// Generated code for Python source for module 'helpers'
// created by Nuitka version 0.5.5.3

// This code is in part copyright 2014 Kay Hayen.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "nuitka/prelude.hpp"

#include "__helpers.hpp"

// The _module_helpers is a Python object pointer of module type.

// Note: For full compatability with CPython, every module variable access
// needs to go through it except for cases where the module cannot possibly
// have changed in the mean time.

PyObject *module_helpers;
PyDictObject *moduledict_helpers;

// The module constants used
extern PyObject *const_int_0;
static PyObject *const_float_0_0;
static PyObject *const_float_1_0;
static PyObject *const_int_neg_1;
extern PyObject *const_int_pos_1;
static PyObject *const_int_pos_2;
static PyObject *const_int_pos_3;
extern PyObject *const_dict_empty;
static PyObject *const_str_plain_a;
static PyObject *const_str_plain_b;
static PyObject *const_str_plain_i;
static PyObject *const_str_plain_n;
static PyObject *const_str_plain_t;
static PyObject *const_str_plain_u;
static PyObject *const_str_plain_v;
extern PyObject *const_tuple_empty;
static PyObject *const_str_plain_dv;
static PyObject *const_str_plain_px;
static PyObject *const_str_plain_py;
static PyObject *const_str_plain_pz;
static PyObject *const_str_plain_qx;
static PyObject *const_str_plain_qy;
static PyObject *const_str_plain_vx;
static PyObject *const_str_plain_vy;
static PyObject *const_str_plain_vz;
static PyObject *const_str_plain_abs;
static PyObject *const_str_plain_dpx;
static PyObject *const_str_plain_dpy;
static PyObject *const_str_plain_dpz;
static PyObject *const_str_plain_dqx;
static PyObject *const_str_plain_dqy;
static PyObject *const_str_plain_dvx;
static PyObject *const_str_plain_dvy;
static PyObject *const_str_plain_dvz;
static PyObject *const_str_plain_pow;
static PyObject *const_str_plain_sign;
static PyObject *const_str_plain_sqrt;
static PyObject *const_str_plain_csqrt;
static PyObject *const_str_plain_floor;
static PyObject *const_str_plain_jumps;
static PyObject *const_str_plain_numpy;
static PyObject *const_str_plain_point;
static PyObject *const_float__minus_1_0;
static PyObject *const_str_plain_append;
static PyObject *const_str_plain_jump_x;
static PyObject *const_str_plain_jump_y;
static PyObject *const_str_plain_jump_z;
static PyObject *const_str_plain_points;
static PyObject *const_float_1e_minus_05;
extern PyObject *const_str_plain___doc__;
static PyObject *const_str_plain_helpers;
extern PyObject *const_str_plain___file__;
static PyObject *const_str_plain_distance;
static PyObject *const_str_plain_crossings;
static PyObject *const_str_plain_intersect;
static PyObject *const_str_plain_jump_mode;
static PyObject *const_str_plain_num_jumps;
static PyObject *const_str_plain_next_point;
static PyObject *const_str_plain_intersect_i;
static PyObject *const_str_plain_intersect_j;
static PyObject *const_str_plain_cross_product;
static PyObject *const_str_plain_crossing_sign;
static PyObject *const_str_plain_current_index;
static PyObject *const_tuple_str_plain_a_tuple;
static PyObject *const_str_plain_already_jumped;
static PyObject *const_str_plain_find_crossings;
static PyObject *const_str_plain_mag_difference;
static PyObject *const_str_plain_segment_lengths;
static PyObject *const_str_plain_comparison_index;
static PyObject *const_str_plain_crossing_direction;
static PyObject *const_str_plain_distance_travelled;
static PyObject *const_str_plain_max_segment_length;
static PyObject *const_str_plain_do_vectors_intersect;
static PyObject *const_str_plain_twice_max_segment_length;
static PyObject *const_tuple_str_plain_a_str_plain_b_tuple;
static PyObject *const_tuple_int_0_float_0_0_float_0_0_tuple;
static PyObject *const_str_digest_3a2dc88a511759902b6efe905e17a85b;
static PyObject *const_str_digest_4524f30e0c6d9d9f108a2d6de2319892;
static PyObject *const_str_digest_65401906cd2c069cbdb366a8b9d824a4;
static PyObject *const_str_digest_b41289c3bcec03f4b67da6ac9f622993;
static PyObject *const_str_digest_e8b4e4b3fe31ecb310658b8a0d0f0c8e;
static PyObject *const_str_digest_f3dba3f4d453eccf489afa145311fe8b;
static PyObject *const_tuple_55a0923f99cdbeabca9978c36dce6b26_tuple;
static PyObject *const_tuple_bd902d90d5066b303273ddc8acbda6e8_tuple;
static PyObject *const_tuple_ef9ce71ce6573efa21bec6d87780d0b0_tuple;
static PyObject *const_tuple_f52371db3c74ce8e796e694ef3cce50d_tuple;
static PyObject *const_tuple_int_0_float__minus_1_0_float__minus_1_0_tuple;
static PyObject *const_tuple_str_plain_px_str_plain_py_str_plain_qx_str_plain_qy_tuple;

static void _initModuleConstants(void)
{
    const_float_0_0 = UNSTREAM_FLOAT( &constant_bin[ 2237 ] );
    const_float_1_0 = UNSTREAM_FLOAT( &constant_bin[ 2245 ] );
    const_int_neg_1 = PyInt_FromLong( -1l );
    const_int_pos_2 = PyInt_FromLong( 2l );
    const_int_pos_3 = PyInt_FromLong( 3l );
    const_str_plain_a = UNSTREAM_CHAR( 97, 1 );
    const_str_plain_b = UNSTREAM_CHAR( 98, 1 );
    const_str_plain_i = UNSTREAM_CHAR( 105, 1 );
    const_str_plain_n = UNSTREAM_CHAR( 110, 1 );
    const_str_plain_t = UNSTREAM_CHAR( 116, 1 );
    const_str_plain_u = UNSTREAM_CHAR( 117, 1 );
    const_str_plain_v = UNSTREAM_CHAR( 118, 1 );
    const_str_plain_dv = UNSTREAM_STRING( &constant_bin[ 86 ], 2, 1 );
    const_str_plain_px = UNSTREAM_STRING( &constant_bin[ 1678 ], 2, 1 );
    const_str_plain_py = UNSTREAM_STRING( &constant_bin[ 1727 ], 2, 1 );
    const_str_plain_pz = UNSTREAM_STRING( &constant_bin[ 2253 ], 2, 1 );
    const_str_plain_qx = UNSTREAM_STRING( &constant_bin[ 1776 ], 2, 1 );
    const_str_plain_qy = UNSTREAM_STRING( &constant_bin[ 1825 ], 2, 1 );
    const_str_plain_vx = UNSTREAM_STRING( &constant_bin[ 491 ], 2, 1 );
    const_str_plain_vy = UNSTREAM_STRING( &constant_bin[ 539 ], 2, 1 );
    const_str_plain_vz = UNSTREAM_STRING( &constant_bin[ 825 ], 2, 1 );
    const_str_plain_abs = UNSTREAM_STRING( &constant_bin[ 2255 ], 3, 1 );
    const_str_plain_dpx = UNSTREAM_STRING( &constant_bin[ 1677 ], 3, 1 );
    const_str_plain_dpy = UNSTREAM_STRING( &constant_bin[ 1726 ], 3, 1 );
    const_str_plain_dpz = UNSTREAM_STRING( &constant_bin[ 2258 ], 3, 1 );
    const_str_plain_dqx = UNSTREAM_STRING( &constant_bin[ 1775 ], 3, 1 );
    const_str_plain_dqy = UNSTREAM_STRING( &constant_bin[ 1824 ], 3, 1 );
    const_str_plain_dvx = UNSTREAM_STRING( &constant_bin[ 587 ], 3, 1 );
    const_str_plain_dvy = UNSTREAM_STRING( &constant_bin[ 636 ], 3, 1 );
    const_str_plain_dvz = UNSTREAM_STRING( &constant_bin[ 930 ], 3, 1 );
    const_str_plain_pow = UNSTREAM_STRING( &constant_bin[ 2261 ], 3, 1 );
    const_str_plain_sign = UNSTREAM_STRING( &constant_bin[ 62 ], 4, 1 );
    const_str_plain_sqrt = UNSTREAM_STRING( &constant_bin[ 295 ], 4, 1 );
    const_str_plain_csqrt = UNSTREAM_STRING( &constant_bin[ 294 ], 5, 1 );
    const_str_plain_floor = UNSTREAM_STRING( &constant_bin[ 1410 ], 5, 1 );
    const_str_plain_jumps = UNSTREAM_STRING( &constant_bin[ 1565 ], 5, 1 );
    const_str_plain_numpy = UNSTREAM_STRING( &constant_bin[ 2264 ], 5, 1 );
    const_str_plain_point = UNSTREAM_STRING( &constant_bin[ 245 ], 5, 1 );
    const_float__minus_1_0 = UNSTREAM_FLOAT( &constant_bin[ 2269 ] );
    const_str_plain_append = UNSTREAM_STRING( &constant_bin[ 2277 ], 6, 1 );
    const_str_plain_jump_x = UNSTREAM_STRING( &constant_bin[ 1078 ], 6, 1 );
    const_str_plain_jump_y = UNSTREAM_STRING( &constant_bin[ 1130 ], 6, 1 );
    const_str_plain_jump_z = UNSTREAM_STRING( &constant_bin[ 740 ], 6, 1 );
    const_str_plain_points = UNSTREAM_STRING( &constant_bin[ 245 ], 6, 1 );
    const_float_1e_minus_05 = UNSTREAM_FLOAT( &constant_bin[ 2283 ] );
    const_str_plain_helpers = UNSTREAM_STRING( &constant_bin[ 2291 ], 7, 1 );
    const_str_plain_distance = UNSTREAM_STRING( &constant_bin[ 1447 ], 8, 1 );
    const_str_plain_crossings = UNSTREAM_STRING( &constant_bin[ 1182 ], 9, 1 );
    const_str_plain_intersect = UNSTREAM_STRING( &constant_bin[ 450 ], 9, 1 );
    const_str_plain_jump_mode = UNSTREAM_STRING( &constant_bin[ 1358 ], 9, 1 );
    const_str_plain_num_jumps = UNSTREAM_STRING( &constant_bin[ 2298 ], 9, 1 );
    const_str_plain_next_point = UNSTREAM_STRING( &constant_bin[ 2307 ], 10, 1 );
    const_str_plain_intersect_i = UNSTREAM_STRING( &constant_bin[ 873 ], 11, 1 );
    const_str_plain_intersect_j = UNSTREAM_STRING( &constant_bin[ 979 ], 11, 1 );
    const_str_plain_cross_product = UNSTREAM_STRING( &constant_bin[ 1033 ], 13, 1 );
    const_str_plain_crossing_sign = UNSTREAM_STRING( &constant_bin[ 2317 ], 13, 1 );
    const_str_plain_current_index = UNSTREAM_STRING( &constant_bin[ 1237 ], 13, 1 );
    const_tuple_str_plain_a_tuple = PyTuple_New( 1 );
    PyTuple_SET_ITEM( const_tuple_str_plain_a_tuple, 0, const_str_plain_a ); Py_INCREF( const_str_plain_a );
    const_str_plain_already_jumped = UNSTREAM_STRING( &constant_bin[ 331 ], 14, 1 );
    const_str_plain_find_crossings = UNSTREAM_STRING( &constant_bin[ 2330 ], 14, 1 );
    const_str_plain_mag_difference = UNSTREAM_STRING( &constant_bin[ 2344 ], 14, 1 );
    const_str_plain_segment_lengths = UNSTREAM_STRING( &constant_bin[ 1616 ], 15, 1 );
    const_str_plain_comparison_index = UNSTREAM_STRING( &constant_bin[ 1296 ], 16, 1 );
    const_str_plain_crossing_direction = UNSTREAM_STRING( &constant_bin[ 2358 ], 18, 1 );
    const_str_plain_distance_travelled = UNSTREAM_STRING( &constant_bin[ 1501 ], 18, 1 );
    const_str_plain_max_segment_length = UNSTREAM_STRING( &constant_bin[ 134 ], 18, 1 );
    const_str_plain_do_vectors_intersect = UNSTREAM_STRING( &constant_bin[ 439 ], 20, 1 );
    const_str_plain_twice_max_segment_length = UNSTREAM_STRING( &constant_bin[ 2376 ], 24, 1 );
    const_tuple_str_plain_a_str_plain_b_tuple = PyTuple_New( 2 );
    PyTuple_SET_ITEM( const_tuple_str_plain_a_str_plain_b_tuple, 0, const_str_plain_a ); Py_INCREF( const_str_plain_a );
    PyTuple_SET_ITEM( const_tuple_str_plain_a_str_plain_b_tuple, 1, const_str_plain_b ); Py_INCREF( const_str_plain_b );
    const_tuple_int_0_float_0_0_float_0_0_tuple = PyTuple_New( 3 );
    PyTuple_SET_ITEM( const_tuple_int_0_float_0_0_float_0_0_tuple, 0, const_int_0 ); Py_INCREF( const_int_0 );
    PyTuple_SET_ITEM( const_tuple_int_0_float_0_0_float_0_0_tuple, 1, const_float_0_0 ); Py_INCREF( const_float_0_0 );
    PyTuple_SET_ITEM( const_tuple_int_0_float_0_0_float_0_0_tuple, 2, const_float_0_0 ); Py_INCREF( const_float_0_0 );
    const_str_digest_3a2dc88a511759902b6efe905e17a85b = UNSTREAM_STRING( &constant_bin[ 2400 ], 184, 0 );
    const_str_digest_4524f30e0c6d9d9f108a2d6de2319892 = UNSTREAM_STRING( &constant_bin[ 2584 ], 838, 0 );
    const_str_digest_65401906cd2c069cbdb366a8b9d824a4 = UNSTREAM_STRING( &constant_bin[ 3422 ], 43, 0 );
    const_str_digest_b41289c3bcec03f4b67da6ac9f622993 = UNSTREAM_STRING( &constant_bin[ 3465 ], 62, 0 );
    const_str_digest_e8b4e4b3fe31ecb310658b8a0d0f0c8e = UNSTREAM_STRING( &constant_bin[ 3527 ], 10, 0 );
    const_str_digest_f3dba3f4d453eccf489afa145311fe8b = UNSTREAM_STRING( &constant_bin[ 3537 ], 43, 0 );
    const_tuple_55a0923f99cdbeabca9978c36dce6b26_tuple = PyTuple_New( 8 );
    PyTuple_SET_ITEM( const_tuple_55a0923f99cdbeabca9978c36dce6b26_tuple, 0, const_str_plain_px ); Py_INCREF( const_str_plain_px );
    PyTuple_SET_ITEM( const_tuple_55a0923f99cdbeabca9978c36dce6b26_tuple, 1, const_str_plain_py ); Py_INCREF( const_str_plain_py );
    PyTuple_SET_ITEM( const_tuple_55a0923f99cdbeabca9978c36dce6b26_tuple, 2, const_str_plain_dpx ); Py_INCREF( const_str_plain_dpx );
    PyTuple_SET_ITEM( const_tuple_55a0923f99cdbeabca9978c36dce6b26_tuple, 3, const_str_plain_dpy ); Py_INCREF( const_str_plain_dpy );
    PyTuple_SET_ITEM( const_tuple_55a0923f99cdbeabca9978c36dce6b26_tuple, 4, const_str_plain_qx ); Py_INCREF( const_str_plain_qx );
    PyTuple_SET_ITEM( const_tuple_55a0923f99cdbeabca9978c36dce6b26_tuple, 5, const_str_plain_qy ); Py_INCREF( const_str_plain_qy );
    PyTuple_SET_ITEM( const_tuple_55a0923f99cdbeabca9978c36dce6b26_tuple, 6, const_str_plain_dqx ); Py_INCREF( const_str_plain_dqx );
    PyTuple_SET_ITEM( const_tuple_55a0923f99cdbeabca9978c36dce6b26_tuple, 7, const_str_plain_dqy ); Py_INCREF( const_str_plain_dqy );
    const_tuple_bd902d90d5066b303273ddc8acbda6e8_tuple = PyMarshal_ReadObjectFromString( (char *)&constant_bin[ 3580 ], 458 );
    const_tuple_ef9ce71ce6573efa21bec6d87780d0b0_tuple = PyTuple_New( 8 );
    PyTuple_SET_ITEM( const_tuple_ef9ce71ce6573efa21bec6d87780d0b0_tuple, 0, const_str_plain_v ); Py_INCREF( const_str_plain_v );
    PyTuple_SET_ITEM( const_tuple_ef9ce71ce6573efa21bec6d87780d0b0_tuple, 1, const_str_plain_dv ); Py_INCREF( const_str_plain_dv );
    PyTuple_SET_ITEM( const_tuple_ef9ce71ce6573efa21bec6d87780d0b0_tuple, 2, const_str_plain_points ); Py_INCREF( const_str_plain_points );
    PyTuple_SET_ITEM( const_tuple_ef9ce71ce6573efa21bec6d87780d0b0_tuple, 3, const_str_plain_segment_lengths ); Py_INCREF( const_str_plain_segment_lengths );
    PyTuple_SET_ITEM( const_tuple_ef9ce71ce6573efa21bec6d87780d0b0_tuple, 4, const_str_plain_current_index ); Py_INCREF( const_str_plain_current_index );
    PyTuple_SET_ITEM( const_tuple_ef9ce71ce6573efa21bec6d87780d0b0_tuple, 5, const_str_plain_comparison_index ); Py_INCREF( const_str_plain_comparison_index );
    PyTuple_SET_ITEM( const_tuple_ef9ce71ce6573efa21bec6d87780d0b0_tuple, 6, const_str_plain_max_segment_length ); Py_INCREF( const_str_plain_max_segment_length );
    PyTuple_SET_ITEM( const_tuple_ef9ce71ce6573efa21bec6d87780d0b0_tuple, 7, const_str_plain_jump_mode ); Py_INCREF( const_str_plain_jump_mode );
    const_tuple_f52371db3c74ce8e796e694ef3cce50d_tuple = PyTuple_New( 10 );
    PyTuple_SET_ITEM( const_tuple_f52371db3c74ce8e796e694ef3cce50d_tuple, 0, const_str_plain_px ); Py_INCREF( const_str_plain_px );
    PyTuple_SET_ITEM( const_tuple_f52371db3c74ce8e796e694ef3cce50d_tuple, 1, const_str_plain_py ); Py_INCREF( const_str_plain_py );
    PyTuple_SET_ITEM( const_tuple_f52371db3c74ce8e796e694ef3cce50d_tuple, 2, const_str_plain_dpx ); Py_INCREF( const_str_plain_dpx );
    PyTuple_SET_ITEM( const_tuple_f52371db3c74ce8e796e694ef3cce50d_tuple, 3, const_str_plain_dpy ); Py_INCREF( const_str_plain_dpy );
    PyTuple_SET_ITEM( const_tuple_f52371db3c74ce8e796e694ef3cce50d_tuple, 4, const_str_plain_qx ); Py_INCREF( const_str_plain_qx );
    PyTuple_SET_ITEM( const_tuple_f52371db3c74ce8e796e694ef3cce50d_tuple, 5, const_str_plain_qy ); Py_INCREF( const_str_plain_qy );
    PyTuple_SET_ITEM( const_tuple_f52371db3c74ce8e796e694ef3cce50d_tuple, 6, const_str_plain_dqx ); Py_INCREF( const_str_plain_dqx );
    PyTuple_SET_ITEM( const_tuple_f52371db3c74ce8e796e694ef3cce50d_tuple, 7, const_str_plain_dqy ); Py_INCREF( const_str_plain_dqy );
    PyTuple_SET_ITEM( const_tuple_f52371db3c74ce8e796e694ef3cce50d_tuple, 8, const_str_plain_t ); Py_INCREF( const_str_plain_t );
    PyTuple_SET_ITEM( const_tuple_f52371db3c74ce8e796e694ef3cce50d_tuple, 9, const_str_plain_u ); Py_INCREF( const_str_plain_u );
    const_tuple_int_0_float__minus_1_0_float__minus_1_0_tuple = PyTuple_New( 3 );
    PyTuple_SET_ITEM( const_tuple_int_0_float__minus_1_0_float__minus_1_0_tuple, 0, const_int_0 ); Py_INCREF( const_int_0 );
    PyTuple_SET_ITEM( const_tuple_int_0_float__minus_1_0_float__minus_1_0_tuple, 1, const_float__minus_1_0 ); Py_INCREF( const_float__minus_1_0 );
    PyTuple_SET_ITEM( const_tuple_int_0_float__minus_1_0_float__minus_1_0_tuple, 2, const_float__minus_1_0 ); Py_INCREF( const_float__minus_1_0 );
    const_tuple_str_plain_px_str_plain_py_str_plain_qx_str_plain_qy_tuple = PyTuple_New( 4 );
    PyTuple_SET_ITEM( const_tuple_str_plain_px_str_plain_py_str_plain_qx_str_plain_qy_tuple, 0, const_str_plain_px ); Py_INCREF( const_str_plain_px );
    PyTuple_SET_ITEM( const_tuple_str_plain_px_str_plain_py_str_plain_qx_str_plain_qy_tuple, 1, const_str_plain_py ); Py_INCREF( const_str_plain_py );
    PyTuple_SET_ITEM( const_tuple_str_plain_px_str_plain_py_str_plain_qx_str_plain_qy_tuple, 2, const_str_plain_qx ); Py_INCREF( const_str_plain_qx );
    PyTuple_SET_ITEM( const_tuple_str_plain_px_str_plain_py_str_plain_qx_str_plain_qy_tuple, 3, const_str_plain_qy ); Py_INCREF( const_str_plain_qy );
}

// The module code objects.
static PyCodeObject *codeobj_641b42b44ac395c476fb706d09df8c7a;
static PyCodeObject *codeobj_38efff3ec8c2b4fd2a37e20a63161fd5;
static PyCodeObject *codeobj_7df00fdabc49ad9593c404a7f3e69637;
static PyCodeObject *codeobj_4f24441a565a494278c4dd4517d38204;
static PyCodeObject *codeobj_ca92ac72b7622f4242ad42d6ae97b5d7;
static PyCodeObject *codeobj_c02aa247315f8c1fd9f538dca91c1671;
static PyCodeObject *codeobj_71cfbdc35b90f7d090c4273357868645;
static PyCodeObject *codeobj_fd2d6789440b495c075e8467bf2253e1;

static void _initModuleCodeObjects(void)
{
    codeobj_641b42b44ac395c476fb706d09df8c7a = MAKE_CODEOBJ( const_str_digest_e8b4e4b3fe31ecb310658b8a0d0f0c8e, const_str_plain_cross_product, 140, const_tuple_str_plain_px_str_plain_py_str_plain_qx_str_plain_qy_tuple, 4, CO_NEWLOCALS | CO_OPTIMIZED | CO_NOFREE );
    codeobj_38efff3ec8c2b4fd2a37e20a63161fd5 = MAKE_CODEOBJ( const_str_digest_e8b4e4b3fe31ecb310658b8a0d0f0c8e, const_str_plain_do_vectors_intersect, 121, const_tuple_55a0923f99cdbeabca9978c36dce6b26_tuple, 8, CO_NEWLOCALS | CO_OPTIMIZED | CO_NOFREE );
    codeobj_7df00fdabc49ad9593c404a7f3e69637 = MAKE_CODEOBJ( const_str_digest_e8b4e4b3fe31ecb310658b8a0d0f0c8e, const_str_plain_do_vectors_intersect, 121, const_tuple_f52371db3c74ce8e796e694ef3cce50d_tuple, 8, CO_NEWLOCALS | CO_OPTIMIZED | CO_NOFREE );
    codeobj_4f24441a565a494278c4dd4517d38204 = MAKE_CODEOBJ( const_str_digest_e8b4e4b3fe31ecb310658b8a0d0f0c8e, const_str_plain_find_crossings, 11, const_tuple_ef9ce71ce6573efa21bec6d87780d0b0_tuple, 8, CO_NEWLOCALS | CO_OPTIMIZED | CO_NOFREE );
    codeobj_ca92ac72b7622f4242ad42d6ae97b5d7 = MAKE_CODEOBJ( const_str_digest_e8b4e4b3fe31ecb310658b8a0d0f0c8e, const_str_plain_find_crossings, 11, const_tuple_bd902d90d5066b303273ddc8acbda6e8_tuple, 8, CO_NEWLOCALS | CO_OPTIMIZED | CO_NOFREE );
    codeobj_c02aa247315f8c1fd9f538dca91c1671 = MAKE_CODEOBJ( const_str_digest_e8b4e4b3fe31ecb310658b8a0d0f0c8e, const_str_plain_helpers, 0, const_tuple_empty, 0, CO_NOFREE );
    codeobj_71cfbdc35b90f7d090c4273357868645 = MAKE_CODEOBJ( const_str_digest_e8b4e4b3fe31ecb310658b8a0d0f0c8e, const_str_plain_mag_difference, 149, const_tuple_str_plain_a_str_plain_b_tuple, 2, CO_NEWLOCALS | CO_OPTIMIZED | CO_NOFREE );
    codeobj_fd2d6789440b495c075e8467bf2253e1 = MAKE_CODEOBJ( const_str_digest_e8b4e4b3fe31ecb310658b8a0d0f0c8e, const_str_plain_sign, 145, const_tuple_str_plain_a_tuple, 1, CO_NEWLOCALS | CO_OPTIMIZED | CO_NOFREE );
}

// The module function declarations.
static PyObject *MAKE_FUNCTION_function_1_find_crossings_of_module_helpers(  );


static PyObject *MAKE_FUNCTION_function_2_do_vectors_intersect_of_module_helpers(  );


static PyObject *MAKE_FUNCTION_function_3_cross_product_of_module_helpers(  );


static PyObject *MAKE_FUNCTION_function_4_sign_of_module_helpers(  );


static PyObject *MAKE_FUNCTION_function_5_mag_difference_of_module_helpers(  );


// The module function definitions.
static PyObject *impl_function_1_find_crossings_of_module_helpers( Nuitka_FunctionObject *self, PyObject *_python_par_v, PyObject *_python_par_dv, PyObject *_python_par_points, PyObject *_python_par_segment_lengths, PyObject *_python_par_current_index, PyObject *_python_par_comparison_index, PyObject *_python_par_max_segment_length, PyObject *_python_par_jump_mode )
{
    // No context is used.

    // Local variable declarations.
    PyObjectLocalVariable par_v; par_v.object = _python_par_v;
    PyObjectLocalVariable par_dv; par_dv.object = _python_par_dv;
    PyObjectLocalVariable par_points; par_points.object = _python_par_points;
    PyObjectLocalVariable par_segment_lengths; par_segment_lengths.object = _python_par_segment_lengths;
    PyObjectLocalVariable par_current_index; par_current_index.object = _python_par_current_index;
    PyObjectLocalVariable par_comparison_index; par_comparison_index.object = _python_par_comparison_index;
    PyObjectLocalVariable par_max_segment_length; par_max_segment_length.object = _python_par_max_segment_length;
    PyObjectLocalVariable par_jump_mode; par_jump_mode.object = _python_par_jump_mode;
    PyObjectLocalVariable var_crossings;
    PyObjectLocalVariable var_vx;
    PyObjectLocalVariable var_vy;
    PyObjectLocalVariable var_vz;
    PyObjectLocalVariable var_dvx;
    PyObjectLocalVariable var_dvy;
    PyObjectLocalVariable var_dvz;
    PyObjectLocalVariable var_twice_max_segment_length;
    PyObjectLocalVariable var_i;
    PyObjectLocalVariable var_already_jumped;
    PyObjectLocalVariable var_point;
    PyObjectLocalVariable var_distance;
    PyObjectLocalVariable var_next_point;
    PyObjectLocalVariable var_jump_x;
    PyObjectLocalVariable var_jump_y;
    PyObjectLocalVariable var_jump_z;
    PyObjectLocalVariable var_intersect;
    PyObjectLocalVariable var_intersect_i;
    PyObjectLocalVariable var_intersect_j;
    PyObjectLocalVariable var_pz;
    PyObjectLocalVariable var_dpz;
    PyObjectLocalVariable var_crossing_sign;
    PyObjectLocalVariable var_crossing_direction;
    PyObjectLocalVariable var_num_jumps;
    PyObjectLocalVariable var_distance_travelled;
    PyObjectLocalVariable var_jumps;
    PyObjectTempVariable tmp_or_1__value_1;
    PyObjectTempVariable tmp_tuple_unpack_1__source_iter;
    PyObjectTempVariable tmp_tuple_unpack_1__element_1;
    PyObjectTempVariable tmp_tuple_unpack_1__element_2;
    PyObjectTempVariable tmp_tuple_unpack_1__element_3;
    PyObjectTempVariable tmp_inplace_assign_1__inplace_start;
    PyObjectTempVariable tmp_inplace_assign_1__inplace_end;
    PyObjectTempVariable tmp_inplace_assign_2__inplace_start;
    PyObjectTempVariable tmp_inplace_assign_2__inplace_end;
    PyObjectTempVariable tmp_inplace_assign_3__inplace_start;
    PyObjectTempVariable tmp_inplace_assign_3__inplace_end;
    PyObjectTempVariable tmp_inplace_assign_4__inplace_start;
    PyObjectTempVariable tmp_inplace_assign_4__inplace_end;
    PyObjectTempVariable tmp_inplace_assign_5__inplace_start;
    PyObjectTempVariable tmp_inplace_assign_5__inplace_end;
    PyObjectTempVariable tmp_inplace_assign_6__inplace_start;
    PyObjectTempVariable tmp_inplace_assign_6__inplace_end;
    PyObjectTempVariable tmp_and_1__value_1;
    PyObjectTempVariable tmp_inplace_assign_7__inplace_start;
    PyObjectTempVariable tmp_inplace_assign_7__inplace_end;
    PyObject *exception_type = NULL, *exception_value = NULL;
    PyTracebackObject *exception_tb = NULL;
    PyObject *exception_keeper_type_1;
    PyObject *exception_keeper_value_1;
    PyTracebackObject *exception_keeper_tb_1;
    PyObject *exception_keeper_type_2;
    PyObject *exception_keeper_value_2;
    PyTracebackObject *exception_keeper_tb_2;
    PyObject *exception_keeper_type_3;
    PyObject *exception_keeper_value_3;
    PyTracebackObject *exception_keeper_tb_3;
    PyObject *exception_keeper_type_4;
    PyObject *exception_keeper_value_4;
    PyTracebackObject *exception_keeper_tb_4;
    PyObject *exception_keeper_type_5;
    PyObject *exception_keeper_value_5;
    PyTracebackObject *exception_keeper_tb_5;
    PyObject *exception_keeper_type_6;
    PyObject *exception_keeper_value_6;
    PyTracebackObject *exception_keeper_tb_6;
    PyObject *exception_keeper_type_7;
    PyObject *exception_keeper_value_7;
    PyTracebackObject *exception_keeper_tb_7;
    PyObject *exception_keeper_type_8;
    PyObject *exception_keeper_value_8;
    PyTracebackObject *exception_keeper_tb_8;
    PyObject *exception_keeper_type_9;
    PyObject *exception_keeper_value_9;
    PyTracebackObject *exception_keeper_tb_9;
    PyObject *exception_keeper_type_10;
    PyObject *exception_keeper_value_10;
    PyTracebackObject *exception_keeper_tb_10;
    PyObject *exception_keeper_type_11;
    PyObject *exception_keeper_value_11;
    PyTracebackObject *exception_keeper_tb_11;
    PyObject *exception_keeper_type_12;
    PyObject *exception_keeper_value_12;
    PyTracebackObject *exception_keeper_tb_12;
    PyObject *exception_keeper_type_13;
    PyObject *exception_keeper_value_13;
    PyTracebackObject *exception_keeper_tb_13;
    PyObject *exception_keeper_type_14;
    PyObject *exception_keeper_value_14;
    PyTracebackObject *exception_keeper_tb_14;
    PyObject *tmp_assign_source_1;
    PyObject *tmp_assign_source_2;
    PyObject *tmp_assign_source_3;
    PyObject *tmp_assign_source_4;
    PyObject *tmp_assign_source_5;
    PyObject *tmp_assign_source_6;
    PyObject *tmp_assign_source_7;
    PyObject *tmp_assign_source_8;
    PyObject *tmp_assign_source_9;
    PyObject *tmp_assign_source_10;
    PyObject *tmp_assign_source_11;
    PyObject *tmp_assign_source_12;
    PyObject *tmp_assign_source_13;
    PyObject *tmp_assign_source_14;
    PyObject *tmp_assign_source_15;
    PyObject *tmp_assign_source_16;
    PyObject *tmp_assign_source_17;
    PyObject *tmp_assign_source_18;
    PyObject *tmp_assign_source_19;
    PyObject *tmp_assign_source_20;
    PyObject *tmp_assign_source_21;
    PyObject *tmp_assign_source_22;
    PyObject *tmp_assign_source_23;
    PyObject *tmp_assign_source_24;
    PyObject *tmp_assign_source_25;
    PyObject *tmp_assign_source_26;
    PyObject *tmp_assign_source_27;
    PyObject *tmp_assign_source_28;
    PyObject *tmp_assign_source_29;
    PyObject *tmp_assign_source_30;
    PyObject *tmp_assign_source_31;
    PyObject *tmp_assign_source_32;
    PyObject *tmp_assign_source_33;
    PyObject *tmp_assign_source_34;
    PyObject *tmp_assign_source_35;
    PyObject *tmp_assign_source_36;
    PyObject *tmp_assign_source_37;
    PyObject *tmp_assign_source_38;
    PyObject *tmp_assign_source_39;
    PyObject *tmp_assign_source_40;
    PyObject *tmp_assign_source_41;
    PyObject *tmp_assign_source_42;
    PyObject *tmp_assign_source_43;
    PyObject *tmp_assign_source_44;
    PyObject *tmp_assign_source_45;
    PyObject *tmp_assign_source_46;
    PyObject *tmp_assign_source_47;
    PyObject *tmp_assign_source_48;
    PyObject *tmp_assign_source_49;
    PyObject *tmp_assign_source_50;
    PyObject *tmp_assign_source_51;
    PyObject *tmp_assign_source_52;
    PyObject *tmp_assign_source_53;
    PyObject *tmp_assign_source_54;
    PyObject *tmp_assign_source_55;
    PyObject *tmp_assign_source_56;
    PyObject *tmp_assign_source_57;
    PyObject *tmp_assign_source_58;
    PyObject *tmp_binop_left_1;
    PyObject *tmp_binop_left_2;
    PyObject *tmp_binop_left_3;
    PyObject *tmp_binop_left_4;
    PyObject *tmp_binop_left_5;
    PyObject *tmp_binop_left_6;
    PyObject *tmp_binop_left_7;
    PyObject *tmp_binop_left_8;
    PyObject *tmp_binop_left_9;
    PyObject *tmp_binop_left_10;
    PyObject *tmp_binop_left_11;
    PyObject *tmp_binop_left_12;
    PyObject *tmp_binop_left_13;
    PyObject *tmp_binop_left_14;
    PyObject *tmp_binop_left_15;
    PyObject *tmp_binop_left_16;
    PyObject *tmp_binop_left_17;
    PyObject *tmp_binop_left_18;
    PyObject *tmp_binop_left_19;
    PyObject *tmp_binop_left_20;
    PyObject *tmp_binop_left_21;
    PyObject *tmp_binop_left_22;
    PyObject *tmp_binop_left_23;
    PyObject *tmp_binop_left_24;
    PyObject *tmp_binop_left_25;
    PyObject *tmp_binop_left_26;
    PyObject *tmp_binop_left_27;
    PyObject *tmp_binop_left_28;
    PyObject *tmp_binop_left_29;
    PyObject *tmp_binop_left_30;
    PyObject *tmp_binop_left_31;
    PyObject *tmp_binop_left_32;
    PyObject *tmp_binop_left_33;
    PyObject *tmp_binop_right_1;
    PyObject *tmp_binop_right_2;
    PyObject *tmp_binop_right_3;
    PyObject *tmp_binop_right_4;
    PyObject *tmp_binop_right_5;
    PyObject *tmp_binop_right_6;
    PyObject *tmp_binop_right_7;
    PyObject *tmp_binop_right_8;
    PyObject *tmp_binop_right_9;
    PyObject *tmp_binop_right_10;
    PyObject *tmp_binop_right_11;
    PyObject *tmp_binop_right_12;
    PyObject *tmp_binop_right_13;
    PyObject *tmp_binop_right_14;
    PyObject *tmp_binop_right_15;
    PyObject *tmp_binop_right_16;
    PyObject *tmp_binop_right_17;
    PyObject *tmp_binop_right_18;
    PyObject *tmp_binop_right_19;
    PyObject *tmp_binop_right_20;
    PyObject *tmp_binop_right_21;
    PyObject *tmp_binop_right_22;
    PyObject *tmp_binop_right_23;
    PyObject *tmp_binop_right_24;
    PyObject *tmp_binop_right_25;
    PyObject *tmp_binop_right_26;
    PyObject *tmp_binop_right_27;
    PyObject *tmp_binop_right_28;
    PyObject *tmp_binop_right_29;
    PyObject *tmp_binop_right_30;
    PyObject *tmp_binop_right_31;
    PyObject *tmp_binop_right_32;
    PyObject *tmp_binop_right_33;
    bool tmp_break_1;
    PyObject *tmp_call_arg_element_1;
    PyObject *tmp_call_arg_element_2;
    PyObject *tmp_call_arg_element_3;
    PyObject *tmp_call_arg_element_4;
    PyObject *tmp_call_arg_element_5;
    PyObject *tmp_call_arg_element_6;
    PyObject *tmp_call_arg_element_7;
    PyObject *tmp_call_arg_element_8;
    PyObject *tmp_call_arg_element_9;
    PyObject *tmp_call_arg_element_10;
    PyObject *tmp_call_arg_element_11;
    PyObject *tmp_call_arg_element_12;
    PyObject *tmp_call_arg_element_13;
    PyObject *tmp_call_arg_element_14;
    PyObject *tmp_call_arg_element_15;
    PyObject *tmp_call_arg_element_16;
    PyObject *tmp_call_arg_element_17;
    PyObject *tmp_call_arg_element_18;
    PyObject *tmp_call_arg_element_19;
    PyObject *tmp_call_arg_element_20;
    PyObject *tmp_call_arg_element_21;
    PyObject *tmp_call_arg_element_22;
    PyObject *tmp_called_1;
    PyObject *tmp_called_2;
    PyObject *tmp_called_3;
    PyObject *tmp_called_4;
    PyObject *tmp_called_5;
    PyObject *tmp_called_6;
    PyObject *tmp_called_7;
    PyObject *tmp_called_8;
    PyObject *tmp_called_9;
    PyObject *tmp_called_10;
    int tmp_cmp_Eq_1;
    int tmp_cmp_Eq_2;
    int tmp_cmp_Gt_1;
    int tmp_cmp_Lt_1;
    int tmp_cmp_Lt_2;
    PyObject *tmp_compare_left_1;
    PyObject *tmp_compare_left_2;
    PyObject *tmp_compare_left_3;
    PyObject *tmp_compare_left_4;
    PyObject *tmp_compare_left_5;
    PyObject *tmp_compare_left_6;
    PyObject *tmp_compare_left_7;
    PyObject *tmp_compare_left_8;
    PyObject *tmp_compare_left_9;
    PyObject *tmp_compare_left_10;
    PyObject *tmp_compare_left_11;
    PyObject *tmp_compare_left_12;
    PyObject *tmp_compare_right_1;
    PyObject *tmp_compare_right_2;
    PyObject *tmp_compare_right_3;
    PyObject *tmp_compare_right_4;
    PyObject *tmp_compare_right_5;
    PyObject *tmp_compare_right_6;
    PyObject *tmp_compare_right_7;
    PyObject *tmp_compare_right_8;
    PyObject *tmp_compare_right_9;
    PyObject *tmp_compare_right_10;
    PyObject *tmp_compare_right_11;
    PyObject *tmp_compare_right_12;
    PyObject *tmp_compexpr_left_1;
    PyObject *tmp_compexpr_left_2;
    PyObject *tmp_compexpr_left_3;
    PyObject *tmp_compexpr_right_1;
    PyObject *tmp_compexpr_right_2;
    PyObject *tmp_compexpr_right_3;
    int tmp_cond_truth_1;
    int tmp_cond_truth_2;
    int tmp_cond_truth_3;
    int tmp_cond_truth_4;
    int tmp_cond_truth_5;
    PyObject *tmp_cond_value_1;
    PyObject *tmp_cond_value_2;
    PyObject *tmp_cond_value_3;
    PyObject *tmp_cond_value_4;
    PyObject *tmp_cond_value_5;
    PyObject *tmp_frame_locals;
    bool tmp_isnot_1;
    bool tmp_isnot_2;
    bool tmp_isnot_3;
    bool tmp_isnot_4;
    bool tmp_isnot_5;
    bool tmp_isnot_6;
    bool tmp_isnot_7;
    PyObject *tmp_iter_arg_1;
    PyObject *tmp_iterator_attempt_1;
    PyObject *tmp_iterator_name_1;
    PyObject *tmp_len_arg_1;
    PyObject *tmp_len_arg_2;
    PyObject *tmp_list_element_1;
    PyObject *tmp_list_element_2;
    bool tmp_result;
    PyObject *tmp_return_value;
    PyObject *tmp_source_name_1;
    PyObject *tmp_source_name_2;
    PyObject *tmp_subscr_subscript_1;
    PyObject *tmp_subscr_subscript_2;
    PyObject *tmp_subscr_subscript_3;
    PyObject *tmp_subscr_subscript_4;
    PyObject *tmp_subscr_subscript_5;
    PyObject *tmp_subscr_subscript_6;
    PyObject *tmp_subscr_subscript_7;
    PyObject *tmp_subscr_subscript_8;
    PyObject *tmp_subscr_subscript_9;
    PyObject *tmp_subscr_subscript_10;
    PyObject *tmp_subscr_subscript_11;
    PyObject *tmp_subscr_subscript_12;
    PyObject *tmp_subscr_subscript_13;
    PyObject *tmp_subscr_subscript_14;
    PyObject *tmp_subscr_subscript_15;
    PyObject *tmp_subscr_subscript_16;
    PyObject *tmp_subscr_subscript_17;
    PyObject *tmp_subscr_subscript_18;
    PyObject *tmp_subscr_subscript_19;
    PyObject *tmp_subscr_subscript_20;
    PyObject *tmp_subscr_target_1;
    PyObject *tmp_subscr_target_2;
    PyObject *tmp_subscr_target_3;
    PyObject *tmp_subscr_target_4;
    PyObject *tmp_subscr_target_5;
    PyObject *tmp_subscr_target_6;
    PyObject *tmp_subscr_target_7;
    PyObject *tmp_subscr_target_8;
    PyObject *tmp_subscr_target_9;
    PyObject *tmp_subscr_target_10;
    PyObject *tmp_subscr_target_11;
    PyObject *tmp_subscr_target_12;
    PyObject *tmp_subscr_target_13;
    PyObject *tmp_subscr_target_14;
    PyObject *tmp_subscr_target_15;
    PyObject *tmp_subscr_target_16;
    PyObject *tmp_subscr_target_17;
    PyObject *tmp_subscr_target_18;
    PyObject *tmp_subscr_target_19;
    PyObject *tmp_subscr_target_20;
    int tmp_tried_lineno_1;
    int tmp_tried_lineno_2;
    int tmp_tried_lineno_3;
    int tmp_tried_lineno_4;
    int tmp_tried_lineno_5;
    int tmp_tried_lineno_6;
    int tmp_tried_lineno_7;
    int tmp_tried_lineno_8;
    int tmp_tried_lineno_9;
    int tmp_tried_lineno_10;
    PyObject *tmp_unpack_1;
    PyObject *tmp_unpack_2;
    PyObject *tmp_unpack_3;
    NUITKA_MAY_BE_UNUSED PyObject *tmp_unused;
    tmp_return_value = NULL;

    // Actual function code.
    tmp_assign_source_1 = PyList_New( 0 );
    assert( var_crossings.object == NULL );
    var_crossings.object = tmp_assign_source_1;

    static PyFrameObject *cache_frame_function = NULL;
    MAKE_OR_REUSE_FRAME( cache_frame_function, codeobj_4f24441a565a494278c4dd4517d38204, module_helpers );
    PyFrameObject *frame_function = cache_frame_function;

    // Push the new frame as the currently active one.
    pushFrameStack( frame_function );

    // Mark the frame object as in use, ref count 1 will be up for reuse.
    Py_INCREF( frame_function );
    assert( Py_REFCNT( frame_function ) == 2 ); // Frame stack

#if PYTHON_VERSION >= 340
    frame_function->f_executing += 1;
#endif

    // Framed code:
    tmp_subscr_target_1 = par_v.object;

    if ( tmp_subscr_target_1 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 23 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 46;
        goto frame_exception_exit_1;
    }

    tmp_subscr_subscript_1 = const_int_0;
    tmp_assign_source_2 = LOOKUP_SUBSCRIPT( tmp_subscr_target_1, tmp_subscr_subscript_1 );
    if ( tmp_assign_source_2 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 46;
        goto frame_exception_exit_1;
    }
    assert( var_vx.object == NULL );
    var_vx.object = tmp_assign_source_2;

    tmp_subscr_target_2 = par_v.object;

    if ( tmp_subscr_target_2 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 23 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 47;
        goto frame_exception_exit_1;
    }

    tmp_subscr_subscript_2 = const_int_pos_1;
    tmp_assign_source_3 = LOOKUP_SUBSCRIPT( tmp_subscr_target_2, tmp_subscr_subscript_2 );
    if ( tmp_assign_source_3 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 47;
        goto frame_exception_exit_1;
    }
    assert( var_vy.object == NULL );
    var_vy.object = tmp_assign_source_3;

    tmp_subscr_target_3 = par_v.object;

    if ( tmp_subscr_target_3 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 23 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 48;
        goto frame_exception_exit_1;
    }

    tmp_subscr_subscript_3 = const_int_pos_2;
    tmp_assign_source_4 = LOOKUP_SUBSCRIPT( tmp_subscr_target_3, tmp_subscr_subscript_3 );
    if ( tmp_assign_source_4 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 48;
        goto frame_exception_exit_1;
    }
    assert( var_vz.object == NULL );
    var_vz.object = tmp_assign_source_4;

    tmp_subscr_target_4 = par_dv.object;

    if ( tmp_subscr_target_4 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 70 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 49;
        goto frame_exception_exit_1;
    }

    tmp_subscr_subscript_4 = const_int_0;
    tmp_assign_source_5 = LOOKUP_SUBSCRIPT( tmp_subscr_target_4, tmp_subscr_subscript_4 );
    if ( tmp_assign_source_5 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 49;
        goto frame_exception_exit_1;
    }
    assert( var_dvx.object == NULL );
    var_dvx.object = tmp_assign_source_5;

    tmp_subscr_target_5 = par_dv.object;

    if ( tmp_subscr_target_5 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 70 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 50;
        goto frame_exception_exit_1;
    }

    tmp_subscr_subscript_5 = const_int_pos_1;
    tmp_assign_source_6 = LOOKUP_SUBSCRIPT( tmp_subscr_target_5, tmp_subscr_subscript_5 );
    if ( tmp_assign_source_6 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 50;
        goto frame_exception_exit_1;
    }
    assert( var_dvy.object == NULL );
    var_dvy.object = tmp_assign_source_6;

    tmp_subscr_target_6 = par_dv.object;

    if ( tmp_subscr_target_6 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 70 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 51;
        goto frame_exception_exit_1;
    }

    tmp_subscr_subscript_6 = const_int_pos_2;
    tmp_assign_source_7 = LOOKUP_SUBSCRIPT( tmp_subscr_target_6, tmp_subscr_subscript_6 );
    if ( tmp_assign_source_7 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 51;
        goto frame_exception_exit_1;
    }
    assert( var_dvz.object == NULL );
    var_dvz.object = tmp_assign_source_7;

    tmp_binop_left_1 = const_int_pos_2;
    tmp_binop_right_1 = par_max_segment_length.object;

    if ( tmp_binop_right_1 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 118 ], 64, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 53;
        goto frame_exception_exit_1;
    }

    tmp_assign_source_8 = BINARY_OPERATION_MUL( tmp_binop_left_1, tmp_binop_right_1 );
    if ( tmp_assign_source_8 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 53;
        goto frame_exception_exit_1;
    }
    assert( var_twice_max_segment_length.object == NULL );
    var_twice_max_segment_length.object = tmp_assign_source_8;

    tmp_assign_source_9 = const_int_0;
    assert( var_i.object == NULL );
    var_i.object = INCREASE_REFCOUNT( tmp_assign_source_9 );

    tmp_assign_source_10 = const_int_0;
    assert( var_already_jumped.object == NULL );
    var_already_jumped.object = INCREASE_REFCOUNT( tmp_assign_source_10 );

    loop_start_1:;
    tmp_compare_left_1 = var_i.object;

    if ( tmp_compare_left_1 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 182 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 58;
        goto frame_exception_exit_1;
    }

    tmp_len_arg_1 = par_points.object;

    if ( tmp_len_arg_1 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 229 ], 52, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 58;
        goto frame_exception_exit_1;
    }

    tmp_binop_left_2 = BUILTIN_LEN( tmp_len_arg_1 );
    if ( tmp_binop_left_2 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 58;
        goto frame_exception_exit_1;
    }
    tmp_binop_right_2 = const_int_pos_1;
    tmp_compare_right_1 = BINARY_OPERATION_SUB( tmp_binop_left_2, tmp_binop_right_2 );
    Py_DECREF( tmp_binop_left_2 );
    if ( tmp_compare_right_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 58;
        goto frame_exception_exit_1;
    }
    tmp_cmp_Lt_1 = RICH_COMPARE_BOOL_LT( tmp_compare_left_1, tmp_compare_right_1 );
    if ( tmp_cmp_Lt_1 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_compare_right_1 );

        frame_function->f_lineno = 58;
        goto frame_exception_exit_1;
    }
    Py_DECREF( tmp_compare_right_1 );
    if (tmp_cmp_Lt_1 == 1)
    {
        goto branch_no_1;
    }
    else
    {
        goto branch_yes_1;
    }
    branch_yes_1:;
    goto loop_end_1;
    branch_no_1:;
    tmp_subscr_target_7 = par_points.object;

    if ( tmp_subscr_target_7 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 229 ], 52, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 59;
        goto frame_exception_exit_1;
    }

    tmp_subscr_subscript_7 = var_i.object;

    if ( tmp_subscr_subscript_7 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 182 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 59;
        goto frame_exception_exit_1;
    }

    tmp_assign_source_11 = LOOKUP_SUBSCRIPT( tmp_subscr_target_7, tmp_subscr_subscript_7 );
    if ( tmp_assign_source_11 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 59;
        goto frame_exception_exit_1;
    }
    if (var_point.object == NULL)
    {
        var_point.object = tmp_assign_source_11;
    }
    else
    {
        PyObject *old = var_point.object;
        var_point.object = tmp_assign_source_11;
        Py_DECREF( old );
    }
    tmp_called_1 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_csqrt );

    if (unlikely( tmp_called_1 == NULL ))
    {
        tmp_called_1 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_csqrt );
    }

    if ( tmp_called_1 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 281 ], 34, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 60;
        goto frame_exception_exit_1;
    }

    tmp_called_2 = LOOKUP_BUILTIN( const_str_plain_pow );
    if ( tmp_called_2 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 60;
        goto frame_exception_exit_1;
    }
    tmp_binop_left_4 = var_vx.object;

    tmp_subscr_target_8 = var_point.object;

    tmp_subscr_subscript_8 = const_int_0;
    tmp_binop_right_4 = LOOKUP_SUBSCRIPT( tmp_subscr_target_8, tmp_subscr_subscript_8 );
    if ( tmp_binop_right_4 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 60;
        goto frame_exception_exit_1;
    }
    tmp_call_arg_element_2 = BINARY_OPERATION_SUB( tmp_binop_left_4, tmp_binop_right_4 );
    Py_DECREF( tmp_binop_right_4 );
    if ( tmp_call_arg_element_2 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 60;
        goto frame_exception_exit_1;
    }
    tmp_call_arg_element_3 = const_int_pos_2;
    frame_function->f_lineno = 60;
    tmp_binop_left_3 = CALL_FUNCTION_WITH_ARGS2( tmp_called_2, tmp_call_arg_element_2, tmp_call_arg_element_3 );
    Py_DECREF( tmp_call_arg_element_2 );
    if ( tmp_binop_left_3 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 60;
        goto frame_exception_exit_1;
    }
    tmp_called_3 = LOOKUP_BUILTIN( const_str_plain_pow );
    if ( tmp_called_3 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_3 );

        frame_function->f_lineno = 60;
        goto frame_exception_exit_1;
    }
    tmp_binop_left_5 = var_vy.object;

    tmp_subscr_target_9 = var_point.object;

    tmp_subscr_subscript_9 = const_int_pos_1;
    tmp_binop_right_5 = LOOKUP_SUBSCRIPT( tmp_subscr_target_9, tmp_subscr_subscript_9 );
    if ( tmp_binop_right_5 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_3 );

        frame_function->f_lineno = 60;
        goto frame_exception_exit_1;
    }
    tmp_call_arg_element_4 = BINARY_OPERATION_SUB( tmp_binop_left_5, tmp_binop_right_5 );
    Py_DECREF( tmp_binop_right_5 );
    if ( tmp_call_arg_element_4 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_3 );

        frame_function->f_lineno = 60;
        goto frame_exception_exit_1;
    }
    tmp_call_arg_element_5 = const_int_pos_2;
    frame_function->f_lineno = 60;
    tmp_binop_right_3 = CALL_FUNCTION_WITH_ARGS2( tmp_called_3, tmp_call_arg_element_4, tmp_call_arg_element_5 );
    Py_DECREF( tmp_call_arg_element_4 );
    if ( tmp_binop_right_3 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_3 );

        frame_function->f_lineno = 60;
        goto frame_exception_exit_1;
    }
    tmp_call_arg_element_1 = BINARY_OPERATION_ADD( tmp_binop_left_3, tmp_binop_right_3 );
    Py_DECREF( tmp_binop_left_3 );
    Py_DECREF( tmp_binop_right_3 );
    if ( tmp_call_arg_element_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 60;
        goto frame_exception_exit_1;
    }
    frame_function->f_lineno = 60;
    tmp_assign_source_12 = CALL_FUNCTION_WITH_ARGS1( tmp_called_1, tmp_call_arg_element_1 );
    Py_DECREF( tmp_call_arg_element_1 );
    if ( tmp_assign_source_12 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 60;
        goto frame_exception_exit_1;
    }
    if (var_distance.object == NULL)
    {
        var_distance.object = tmp_assign_source_12;
    }
    else
    {
        PyObject *old = var_distance.object;
        var_distance.object = tmp_assign_source_12;
        Py_DECREF( old );
    }
    // Tried code
    tmp_cond_value_1 = NULL;
    // Tried code
    tmp_compexpr_left_1 = var_distance.object;

    tmp_compexpr_right_1 = var_twice_max_segment_length.object;

    tmp_assign_source_13 = RICH_COMPARE_LT( tmp_compexpr_left_1, tmp_compexpr_right_1 );
    if ( tmp_assign_source_13 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 61;
        goto try_finally_handler_2;
    }
    if (tmp_or_1__value_1.object == NULL)
    {
        tmp_or_1__value_1.object = tmp_assign_source_13;
    }
    else
    {
        PyObject *old = tmp_or_1__value_1.object;
        tmp_or_1__value_1.object = tmp_assign_source_13;
        Py_DECREF( old );
    }
    tmp_cond_value_2 = tmp_or_1__value_1.object;

    tmp_cond_truth_2 = CHECK_IF_TRUE( tmp_cond_value_2 );
    if ( tmp_cond_truth_2 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 61;
        goto try_finally_handler_2;
    }
    if (tmp_cond_truth_2 == 1)
    {
        goto condexpr_true_1;
    }
    else
    {
        goto condexpr_false_1;
    }
    condexpr_true_1:;
    tmp_cond_value_1 = tmp_or_1__value_1.object;

    goto condexpr_end_1;
    condexpr_false_1:;
    tmp_cond_value_1 = NULL;
    // Tried code
    tmp_result = tmp_or_1__value_1.object != NULL;
    if ( tmp_result == true )
    {
        Py_DECREF( tmp_or_1__value_1.object );
        tmp_or_1__value_1.object = NULL;
    }

    assert( tmp_result != false );
    tmp_cond_value_1 = var_already_jumped.object;

    if ( tmp_cond_value_1 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 315 ], 60, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 61;
        goto try_finally_handler_3;
    }

    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_3:;
    exception_keeper_type_1 = exception_type;
    exception_keeper_value_1 = exception_value;
    exception_keeper_tb_1 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_1 != NULL )
    {
        exception_type = exception_keeper_type_1;
        exception_value = exception_keeper_value_1;
        exception_tb = exception_keeper_tb_1;

        goto try_finally_handler_2;
    }

    goto finally_end_1;
    finally_end_1:;
    condexpr_end_1:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_2:;
    exception_keeper_type_2 = exception_type;
    exception_keeper_value_2 = exception_value;
    exception_keeper_tb_2 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_2 != NULL )
    {
        exception_type = exception_keeper_type_2;
        exception_value = exception_keeper_value_2;
        exception_tb = exception_keeper_tb_2;

        goto try_finally_handler_1;
    }

    goto finally_end_2;
    finally_end_2:;
    tmp_cond_truth_1 = CHECK_IF_TRUE( tmp_cond_value_1 );
    if ( tmp_cond_truth_1 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 61;
        goto try_finally_handler_1;
    }
    if (tmp_cond_truth_1 == 1)
    {
        goto branch_yes_2;
    }
    else
    {
        goto branch_no_2;
    }
    branch_yes_2:;
    tmp_assign_source_14 = const_int_0;
    if (var_already_jumped.object == NULL)
    {
        var_already_jumped.object = INCREASE_REFCOUNT( tmp_assign_source_14 );
    }
    else
    {
        PyObject *old = var_already_jumped.object;
        var_already_jumped.object = INCREASE_REFCOUNT( tmp_assign_source_14 );
        Py_DECREF( old );
    }
    tmp_subscr_target_10 = par_points.object;

    if ( tmp_subscr_target_10 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 229 ], 52, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 63;
        goto try_finally_handler_1;
    }

    tmp_binop_left_6 = var_i.object;

    if ( tmp_binop_left_6 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 182 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 63;
        goto try_finally_handler_1;
    }

    tmp_binop_right_6 = const_int_pos_1;
    tmp_subscr_subscript_10 = BINARY_OPERATION_ADD( tmp_binop_left_6, tmp_binop_right_6 );
    if ( tmp_subscr_subscript_10 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 63;
        goto try_finally_handler_1;
    }
    tmp_assign_source_15 = LOOKUP_SUBSCRIPT( tmp_subscr_target_10, tmp_subscr_subscript_10 );
    Py_DECREF( tmp_subscr_subscript_10 );
    if ( tmp_assign_source_15 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 63;
        goto try_finally_handler_1;
    }
    if (var_next_point.object == NULL)
    {
        var_next_point.object = tmp_assign_source_15;
    }
    else
    {
        PyObject *old = var_next_point.object;
        var_next_point.object = tmp_assign_source_15;
        Py_DECREF( old );
    }
    tmp_subscr_target_11 = var_next_point.object;

    tmp_subscr_subscript_11 = const_int_0;
    tmp_binop_left_7 = LOOKUP_SUBSCRIPT( tmp_subscr_target_11, tmp_subscr_subscript_11 );
    if ( tmp_binop_left_7 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 64;
        goto try_finally_handler_1;
    }
    tmp_subscr_target_12 = var_point.object;

    if ( tmp_subscr_target_12 == NULL )
    {
        Py_DECREF( tmp_binop_left_7 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 375 ], 51, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 64;
        goto try_finally_handler_1;
    }

    tmp_subscr_subscript_12 = const_int_0;
    tmp_binop_right_7 = LOOKUP_SUBSCRIPT( tmp_subscr_target_12, tmp_subscr_subscript_12 );
    if ( tmp_binop_right_7 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_7 );

        frame_function->f_lineno = 64;
        goto try_finally_handler_1;
    }
    tmp_assign_source_16 = BINARY_OPERATION_SUB( tmp_binop_left_7, tmp_binop_right_7 );
    Py_DECREF( tmp_binop_left_7 );
    Py_DECREF( tmp_binop_right_7 );
    if ( tmp_assign_source_16 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 64;
        goto try_finally_handler_1;
    }
    if (var_jump_x.object == NULL)
    {
        var_jump_x.object = tmp_assign_source_16;
    }
    else
    {
        PyObject *old = var_jump_x.object;
        var_jump_x.object = tmp_assign_source_16;
        Py_DECREF( old );
    }
    tmp_subscr_target_13 = var_next_point.object;

    tmp_subscr_subscript_13 = const_int_pos_1;
    tmp_binop_left_8 = LOOKUP_SUBSCRIPT( tmp_subscr_target_13, tmp_subscr_subscript_13 );
    if ( tmp_binop_left_8 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 65;
        goto try_finally_handler_1;
    }
    tmp_subscr_target_14 = var_point.object;

    if ( tmp_subscr_target_14 == NULL )
    {
        Py_DECREF( tmp_binop_left_8 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 375 ], 51, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 65;
        goto try_finally_handler_1;
    }

    tmp_subscr_subscript_14 = const_int_pos_1;
    tmp_binop_right_8 = LOOKUP_SUBSCRIPT( tmp_subscr_target_14, tmp_subscr_subscript_14 );
    if ( tmp_binop_right_8 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_8 );

        frame_function->f_lineno = 65;
        goto try_finally_handler_1;
    }
    tmp_assign_source_17 = BINARY_OPERATION_SUB( tmp_binop_left_8, tmp_binop_right_8 );
    Py_DECREF( tmp_binop_left_8 );
    Py_DECREF( tmp_binop_right_8 );
    if ( tmp_assign_source_17 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 65;
        goto try_finally_handler_1;
    }
    if (var_jump_y.object == NULL)
    {
        var_jump_y.object = tmp_assign_source_17;
    }
    else
    {
        PyObject *old = var_jump_y.object;
        var_jump_y.object = tmp_assign_source_17;
        Py_DECREF( old );
    }
    tmp_subscr_target_15 = var_next_point.object;

    tmp_subscr_subscript_15 = const_int_pos_2;
    tmp_binop_left_9 = LOOKUP_SUBSCRIPT( tmp_subscr_target_15, tmp_subscr_subscript_15 );
    if ( tmp_binop_left_9 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 66;
        goto try_finally_handler_1;
    }
    tmp_subscr_target_16 = var_point.object;

    if ( tmp_subscr_target_16 == NULL )
    {
        Py_DECREF( tmp_binop_left_9 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 375 ], 51, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 66;
        goto try_finally_handler_1;
    }

    tmp_subscr_subscript_16 = const_int_pos_2;
    tmp_binop_right_9 = LOOKUP_SUBSCRIPT( tmp_subscr_target_16, tmp_subscr_subscript_16 );
    if ( tmp_binop_right_9 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_9 );

        frame_function->f_lineno = 66;
        goto try_finally_handler_1;
    }
    tmp_assign_source_18 = BINARY_OPERATION_SUB( tmp_binop_left_9, tmp_binop_right_9 );
    Py_DECREF( tmp_binop_left_9 );
    Py_DECREF( tmp_binop_right_9 );
    if ( tmp_assign_source_18 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 66;
        goto try_finally_handler_1;
    }
    if (var_jump_z.object == NULL)
    {
        var_jump_z.object = tmp_assign_source_18;
    }
    else
    {
        PyObject *old = var_jump_z.object;
        var_jump_z.object = tmp_assign_source_18;
        Py_DECREF( old );
    }
    // Tried code
    tmp_called_4 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_do_vectors_intersect );

    if (unlikely( tmp_called_4 == NULL ))
    {
        tmp_called_4 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_do_vectors_intersect );
    }

    if ( tmp_called_4 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 426 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 69;
        goto try_finally_handler_4;
    }

    tmp_call_arg_element_6 = var_vx.object;

    if ( tmp_call_arg_element_6 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 475 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 70;
        goto try_finally_handler_4;
    }

    tmp_call_arg_element_7 = var_vy.object;

    if ( tmp_call_arg_element_7 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 523 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 70;
        goto try_finally_handler_4;
    }

    tmp_call_arg_element_8 = var_dvx.object;

    if ( tmp_call_arg_element_8 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 571 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 70;
        goto try_finally_handler_4;
    }

    tmp_call_arg_element_9 = var_dvy.object;

    if ( tmp_call_arg_element_9 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 620 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 70;
        goto try_finally_handler_4;
    }

    tmp_subscr_target_17 = var_point.object;

    if ( tmp_subscr_target_17 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 375 ], 51, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 70;
        goto try_finally_handler_4;
    }

    tmp_subscr_subscript_17 = const_int_0;
    tmp_call_arg_element_10 = LOOKUP_SUBSCRIPT( tmp_subscr_target_17, tmp_subscr_subscript_17 );
    if ( tmp_call_arg_element_10 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 70;
        goto try_finally_handler_4;
    }
    tmp_subscr_target_18 = var_point.object;

    if ( tmp_subscr_target_18 == NULL )
    {
        Py_DECREF( tmp_call_arg_element_10 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 375 ], 51, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 70;
        goto try_finally_handler_4;
    }

    tmp_subscr_subscript_18 = const_int_pos_1;
    tmp_call_arg_element_11 = LOOKUP_SUBSCRIPT( tmp_subscr_target_18, tmp_subscr_subscript_18 );
    if ( tmp_call_arg_element_11 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_call_arg_element_10 );

        frame_function->f_lineno = 70;
        goto try_finally_handler_4;
    }
    tmp_call_arg_element_12 = var_jump_x.object;

    tmp_call_arg_element_13 = var_jump_y.object;

    frame_function->f_lineno = 71;
    tmp_iter_arg_1 = CALL_FUNCTION_WITH_ARGS8( tmp_called_4, tmp_call_arg_element_6, tmp_call_arg_element_7, tmp_call_arg_element_8, tmp_call_arg_element_9, tmp_call_arg_element_10, tmp_call_arg_element_11, tmp_call_arg_element_12, tmp_call_arg_element_13 );
    Py_DECREF( tmp_call_arg_element_10 );
    Py_DECREF( tmp_call_arg_element_11 );
    if ( tmp_iter_arg_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 71;
        goto try_finally_handler_4;
    }
    tmp_assign_source_19 = MAKE_ITERATOR( tmp_iter_arg_1 );
    Py_DECREF( tmp_iter_arg_1 );
    if ( tmp_assign_source_19 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 69;
        goto try_finally_handler_4;
    }
    if (tmp_tuple_unpack_1__source_iter.object == NULL)
    {
        tmp_tuple_unpack_1__source_iter.object = tmp_assign_source_19;
    }
    else
    {
        PyObject *old = tmp_tuple_unpack_1__source_iter.object;
        tmp_tuple_unpack_1__source_iter.object = tmp_assign_source_19;
        Py_DECREF( old );
    }
    tmp_unpack_1 = tmp_tuple_unpack_1__source_iter.object;

    tmp_assign_source_20 = UNPACK_PARAMETER_NEXT( tmp_unpack_1, 0 );
    if ( tmp_assign_source_20 == NULL )
    {
        if ( !ERROR_OCCURED() )
        {
            exception_type = INCREASE_REFCOUNT( PyExc_StopIteration );
            exception_value = NULL;
            exception_tb = NULL;
        }
        else
        {
            PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        }


        frame_function->f_lineno = 69;
        goto try_finally_handler_4;
    }
    if (tmp_tuple_unpack_1__element_1.object == NULL)
    {
        tmp_tuple_unpack_1__element_1.object = tmp_assign_source_20;
    }
    else
    {
        PyObject *old = tmp_tuple_unpack_1__element_1.object;
        tmp_tuple_unpack_1__element_1.object = tmp_assign_source_20;
        Py_DECREF( old );
    }
    tmp_unpack_2 = tmp_tuple_unpack_1__source_iter.object;

    tmp_assign_source_21 = UNPACK_PARAMETER_NEXT( tmp_unpack_2, 1 );
    if ( tmp_assign_source_21 == NULL )
    {
        if ( !ERROR_OCCURED() )
        {
            exception_type = INCREASE_REFCOUNT( PyExc_StopIteration );
            exception_value = NULL;
            exception_tb = NULL;
        }
        else
        {
            PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        }


        frame_function->f_lineno = 69;
        goto try_finally_handler_4;
    }
    if (tmp_tuple_unpack_1__element_2.object == NULL)
    {
        tmp_tuple_unpack_1__element_2.object = tmp_assign_source_21;
    }
    else
    {
        PyObject *old = tmp_tuple_unpack_1__element_2.object;
        tmp_tuple_unpack_1__element_2.object = tmp_assign_source_21;
        Py_DECREF( old );
    }
    tmp_unpack_3 = tmp_tuple_unpack_1__source_iter.object;

    tmp_assign_source_22 = UNPACK_PARAMETER_NEXT( tmp_unpack_3, 2 );
    if ( tmp_assign_source_22 == NULL )
    {
        if ( !ERROR_OCCURED() )
        {
            exception_type = INCREASE_REFCOUNT( PyExc_StopIteration );
            exception_value = NULL;
            exception_tb = NULL;
        }
        else
        {
            PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        }


        frame_function->f_lineno = 69;
        goto try_finally_handler_4;
    }
    if (tmp_tuple_unpack_1__element_3.object == NULL)
    {
        tmp_tuple_unpack_1__element_3.object = tmp_assign_source_22;
    }
    else
    {
        PyObject *old = tmp_tuple_unpack_1__element_3.object;
        tmp_tuple_unpack_1__element_3.object = tmp_assign_source_22;
        Py_DECREF( old );
    }
    tmp_iterator_name_1 = tmp_tuple_unpack_1__source_iter.object;

    // Check if iterator has left-over elements.
    assertObject( tmp_iterator_name_1 ); assert( PyIter_Check( tmp_iterator_name_1 ) );

    tmp_iterator_attempt_1 = (*Py_TYPE( tmp_iterator_name_1 )->tp_iternext)( tmp_iterator_name_1 );

    if (likely( tmp_iterator_attempt_1 == NULL ))
    {
        // TODO: Could first fetch, then check, should be faster.
        if ( !ERROR_OCCURED() )
        {
        }
        else if ( PyErr_ExceptionMatches( PyExc_StopIteration ))
        {
            PyErr_Clear();
        }
        else
        {
            PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );

            goto try_finally_handler_4;
        }
    }
    else
    {
        Py_DECREF( tmp_iterator_attempt_1 );

        // TODO: Could avoid PyErr_Format.
#if PYTHON_VERSION < 300
        PyErr_Format( PyExc_ValueError, "too many values to unpack" );
#else
        PyErr_Format( PyExc_ValueError, "too many values to unpack (expected 3)" );
#endif
        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );

        goto try_finally_handler_4;
    }
    tmp_assign_source_23 = tmp_tuple_unpack_1__element_1.object;

    if (var_intersect.object == NULL)
    {
        var_intersect.object = INCREASE_REFCOUNT( tmp_assign_source_23 );
    }
    else
    {
        PyObject *old = var_intersect.object;
        var_intersect.object = INCREASE_REFCOUNT( tmp_assign_source_23 );
        Py_DECREF( old );
    }
    tmp_assign_source_24 = tmp_tuple_unpack_1__element_2.object;

    if (var_intersect_i.object == NULL)
    {
        var_intersect_i.object = INCREASE_REFCOUNT( tmp_assign_source_24 );
    }
    else
    {
        PyObject *old = var_intersect_i.object;
        var_intersect_i.object = INCREASE_REFCOUNT( tmp_assign_source_24 );
        Py_DECREF( old );
    }
    tmp_assign_source_25 = tmp_tuple_unpack_1__element_3.object;

    if (var_intersect_j.object == NULL)
    {
        var_intersect_j.object = INCREASE_REFCOUNT( tmp_assign_source_25 );
    }
    else
    {
        PyObject *old = var_intersect_j.object;
        var_intersect_j.object = INCREASE_REFCOUNT( tmp_assign_source_25 );
        Py_DECREF( old );
    }
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_4:;
    exception_keeper_type_3 = exception_type;
    exception_keeper_value_3 = exception_value;
    exception_keeper_tb_3 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    tmp_tried_lineno_1 = frame_function->f_lineno;
    Py_XDECREF( tmp_tuple_unpack_1__source_iter.object );
    tmp_tuple_unpack_1__source_iter.object = NULL;

    Py_XDECREF( tmp_tuple_unpack_1__element_1.object );
    tmp_tuple_unpack_1__element_1.object = NULL;

    Py_XDECREF( tmp_tuple_unpack_1__element_2.object );
    tmp_tuple_unpack_1__element_2.object = NULL;

    Py_XDECREF( tmp_tuple_unpack_1__element_3.object );
    tmp_tuple_unpack_1__element_3.object = NULL;

    frame_function->f_lineno = tmp_tried_lineno_1;
    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_3 != NULL )
    {
        exception_type = exception_keeper_type_3;
        exception_value = exception_keeper_value_3;
        exception_tb = exception_keeper_tb_3;

        goto try_finally_handler_1;
    }

    goto finally_end_3;
    finally_end_3:;
    tmp_cond_value_3 = var_intersect.object;

    if ( tmp_cond_value_3 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 669 ], 55, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 73;
        goto try_finally_handler_1;
    }

    tmp_cond_truth_3 = CHECK_IF_TRUE( tmp_cond_value_3 );
    if ( tmp_cond_truth_3 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 73;
        goto try_finally_handler_1;
    }
    if (tmp_cond_truth_3 == 1)
    {
        goto branch_yes_3;
    }
    else
    {
        goto branch_no_3;
    }
    branch_yes_3:;
    tmp_subscr_target_19 = var_point.object;

    if ( tmp_subscr_target_19 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 375 ], 51, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 74;
        goto try_finally_handler_1;
    }

    tmp_subscr_subscript_19 = const_int_pos_2;
    tmp_assign_source_26 = LOOKUP_SUBSCRIPT( tmp_subscr_target_19, tmp_subscr_subscript_19 );
    if ( tmp_assign_source_26 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 74;
        goto try_finally_handler_1;
    }
    if (var_pz.object == NULL)
    {
        var_pz.object = tmp_assign_source_26;
    }
    else
    {
        PyObject *old = var_pz.object;
        var_pz.object = tmp_assign_source_26;
        Py_DECREF( old );
    }
    tmp_assign_source_27 = var_jump_z.object;

    if ( tmp_assign_source_27 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 724 ], 52, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 75;
        goto try_finally_handler_1;
    }

    if (var_dpz.object == NULL)
    {
        var_dpz.object = INCREASE_REFCOUNT( tmp_assign_source_27 );
    }
    else
    {
        PyObject *old = var_dpz.object;
        var_dpz.object = INCREASE_REFCOUNT( tmp_assign_source_27 );
        Py_DECREF( old );
    }
    tmp_called_5 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_sign );

    if (unlikely( tmp_called_5 == NULL ))
    {
        tmp_called_5 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_sign );
    }

    if ( tmp_called_5 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 776 ], 33, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 77;
        goto try_finally_handler_1;
    }

    tmp_binop_left_11 = var_vz.object;

    if ( tmp_binop_left_11 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 809 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 77;
        goto try_finally_handler_1;
    }

    tmp_binop_left_12 = var_intersect_i.object;

    if ( tmp_binop_left_12 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 857 ], 57, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 77;
        goto try_finally_handler_1;
    }

    tmp_binop_right_12 = var_dvz.object;

    if ( tmp_binop_right_12 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 914 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 77;
        goto try_finally_handler_1;
    }

    tmp_binop_right_11 = BINARY_OPERATION_MUL( tmp_binop_left_12, tmp_binop_right_12 );
    if ( tmp_binop_right_11 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 77;
        goto try_finally_handler_1;
    }
    tmp_binop_left_10 = BINARY_OPERATION_ADD( tmp_binop_left_11, tmp_binop_right_11 );
    Py_DECREF( tmp_binop_right_11 );
    if ( tmp_binop_left_10 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 77;
        goto try_finally_handler_1;
    }
    tmp_binop_left_13 = var_pz.object;

    tmp_binop_left_14 = var_intersect_j.object;

    if ( tmp_binop_left_14 == NULL )
    {
        Py_DECREF( tmp_binop_left_10 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 963 ], 57, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 78;
        goto try_finally_handler_1;
    }

    tmp_binop_right_14 = var_dpz.object;

    tmp_binop_right_13 = BINARY_OPERATION_MUL( tmp_binop_left_14, tmp_binop_right_14 );
    if ( tmp_binop_right_13 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_10 );

        frame_function->f_lineno = 78;
        goto try_finally_handler_1;
    }
    tmp_binop_right_10 = BINARY_OPERATION_ADD( tmp_binop_left_13, tmp_binop_right_13 );
    Py_DECREF( tmp_binop_right_13 );
    if ( tmp_binop_right_10 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_10 );

        frame_function->f_lineno = 78;
        goto try_finally_handler_1;
    }
    tmp_call_arg_element_14 = BINARY_OPERATION_SUB( tmp_binop_left_10, tmp_binop_right_10 );
    Py_DECREF( tmp_binop_left_10 );
    Py_DECREF( tmp_binop_right_10 );
    if ( tmp_call_arg_element_14 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 77;
        goto try_finally_handler_1;
    }
    frame_function->f_lineno = 77;
    tmp_assign_source_28 = CALL_FUNCTION_WITH_ARGS1( tmp_called_5, tmp_call_arg_element_14 );
    Py_DECREF( tmp_call_arg_element_14 );
    if ( tmp_assign_source_28 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 77;
        goto try_finally_handler_1;
    }
    if (var_crossing_sign.object == NULL)
    {
        var_crossing_sign.object = tmp_assign_source_28;
    }
    else
    {
        PyObject *old = var_crossing_sign.object;
        var_crossing_sign.object = tmp_assign_source_28;
        Py_DECREF( old );
    }
    tmp_called_6 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_sign );

    if (unlikely( tmp_called_6 == NULL ))
    {
        tmp_called_6 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_sign );
    }

    if ( tmp_called_6 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 776 ], 33, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 80;
        goto try_finally_handler_1;
    }

    tmp_called_7 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_cross_product );

    if (unlikely( tmp_called_7 == NULL ))
    {
        tmp_called_7 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_cross_product );
    }

    if ( tmp_called_7 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1020 ], 42, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 80;
        goto try_finally_handler_1;
    }

    tmp_call_arg_element_16 = var_dvx.object;

    if ( tmp_call_arg_element_16 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 571 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 81;
        goto try_finally_handler_1;
    }

    tmp_call_arg_element_17 = var_dvy.object;

    if ( tmp_call_arg_element_17 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 620 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 81;
        goto try_finally_handler_1;
    }

    tmp_call_arg_element_18 = var_jump_x.object;

    if ( tmp_call_arg_element_18 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1062 ], 52, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 81;
        goto try_finally_handler_1;
    }

    tmp_call_arg_element_19 = var_jump_y.object;

    if ( tmp_call_arg_element_19 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1114 ], 52, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 81;
        goto try_finally_handler_1;
    }

    frame_function->f_lineno = 81;
    tmp_call_arg_element_15 = CALL_FUNCTION_WITH_ARGS4( tmp_called_7, tmp_call_arg_element_16, tmp_call_arg_element_17, tmp_call_arg_element_18, tmp_call_arg_element_19 );
    if ( tmp_call_arg_element_15 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 81;
        goto try_finally_handler_1;
    }
    frame_function->f_lineno = 80;
    tmp_assign_source_29 = CALL_FUNCTION_WITH_ARGS1( tmp_called_6, tmp_call_arg_element_15 );
    Py_DECREF( tmp_call_arg_element_15 );
    if ( tmp_assign_source_29 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 80;
        goto try_finally_handler_1;
    }
    if (var_crossing_direction.object == NULL)
    {
        var_crossing_direction.object = tmp_assign_source_29;
    }
    else
    {
        PyObject *old = var_crossing_direction.object;
        var_crossing_direction.object = tmp_assign_source_29;
        Py_DECREF( old );
    }
    tmp_source_name_1 = var_crossings.object;

    if ( tmp_source_name_1 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1166 ], 55, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 83;
        goto try_finally_handler_1;
    }

    tmp_called_8 = LOOKUP_ATTRIBUTE( tmp_source_name_1, const_str_plain_append );
    if ( tmp_called_8 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 83;
        goto try_finally_handler_1;
    }
    tmp_call_arg_element_20 = PyList_New( 4 );
    tmp_binop_left_15 = par_current_index.object;

    if ( tmp_binop_left_15 == NULL )
    {
        Py_DECREF( tmp_called_8 );
        Py_DECREF( tmp_call_arg_element_20 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1221 ], 59, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 83;
        goto try_finally_handler_1;
    }

    tmp_binop_right_15 = var_intersect_i.object;

    if ( tmp_binop_right_15 == NULL )
    {
        Py_DECREF( tmp_called_8 );
        Py_DECREF( tmp_call_arg_element_20 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 857 ], 57, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 83;
        goto try_finally_handler_1;
    }

    tmp_list_element_1 = BINARY_OPERATION_ADD( tmp_binop_left_15, tmp_binop_right_15 );
    if ( tmp_list_element_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_called_8 );
        Py_DECREF( tmp_call_arg_element_20 );

        frame_function->f_lineno = 83;
        goto try_finally_handler_1;
    }
    PyList_SET_ITEM( tmp_call_arg_element_20, 0, tmp_list_element_1 );
    tmp_binop_left_17 = par_comparison_index.object;

    if ( tmp_binop_left_17 == NULL )
    {
        Py_DECREF( tmp_called_8 );
        Py_DECREF( tmp_call_arg_element_20 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1280 ], 62, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 84;
        goto try_finally_handler_1;
    }

    tmp_binop_right_17 = var_intersect_j.object;

    if ( tmp_binop_right_17 == NULL )
    {
        Py_DECREF( tmp_called_8 );
        Py_DECREF( tmp_call_arg_element_20 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 963 ], 57, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 84;
        goto try_finally_handler_1;
    }

    tmp_binop_left_16 = BINARY_OPERATION_ADD( tmp_binop_left_17, tmp_binop_right_17 );
    if ( tmp_binop_left_16 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_called_8 );
        Py_DECREF( tmp_call_arg_element_20 );

        frame_function->f_lineno = 84;
        goto try_finally_handler_1;
    }
    tmp_binop_right_16 = var_i.object;

    if ( tmp_binop_right_16 == NULL )
    {
        Py_DECREF( tmp_called_8 );
        Py_DECREF( tmp_call_arg_element_20 );
        Py_DECREF( tmp_binop_left_16 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 182 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 85;
        goto try_finally_handler_1;
    }

    tmp_list_element_1 = BINARY_OPERATION_ADD( tmp_binop_left_16, tmp_binop_right_16 );
    Py_DECREF( tmp_binop_left_16 );
    if ( tmp_list_element_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_called_8 );
        Py_DECREF( tmp_call_arg_element_20 );

        frame_function->f_lineno = 84;
        goto try_finally_handler_1;
    }
    PyList_SET_ITEM( tmp_call_arg_element_20, 1, tmp_list_element_1 );
    tmp_list_element_1 = var_crossing_sign.object;

    Py_INCREF( tmp_list_element_1 );
    PyList_SET_ITEM( tmp_call_arg_element_20, 2, tmp_list_element_1 );
    tmp_binop_left_18 = var_crossing_sign.object;

    tmp_binop_right_18 = var_crossing_direction.object;

    tmp_list_element_1 = BINARY_OPERATION_MUL( tmp_binop_left_18, tmp_binop_right_18 );
    if ( tmp_list_element_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_called_8 );
        Py_DECREF( tmp_call_arg_element_20 );

        frame_function->f_lineno = 87;
        goto try_finally_handler_1;
    }
    PyList_SET_ITEM( tmp_call_arg_element_20, 3, tmp_list_element_1 );
    frame_function->f_lineno = 87;
    tmp_unused = CALL_FUNCTION_WITH_ARGS1( tmp_called_8, tmp_call_arg_element_20 );
    Py_DECREF( tmp_called_8 );
    Py_DECREF( tmp_call_arg_element_20 );
    if ( tmp_unused == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 87;
        goto try_finally_handler_1;
    }
    Py_DECREF( tmp_unused );
    tmp_source_name_2 = var_crossings.object;

    if ( tmp_source_name_2 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1166 ], 55, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 88;
        goto try_finally_handler_1;
    }

    tmp_called_9 = LOOKUP_ATTRIBUTE( tmp_source_name_2, const_str_plain_append );
    if ( tmp_called_9 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 88;
        goto try_finally_handler_1;
    }
    tmp_call_arg_element_21 = PyList_New( 4 );
    tmp_binop_left_20 = par_comparison_index.object;

    if ( tmp_binop_left_20 == NULL )
    {
        Py_DECREF( tmp_called_9 );
        Py_DECREF( tmp_call_arg_element_21 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1280 ], 62, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 88;
        goto try_finally_handler_1;
    }

    tmp_binop_right_20 = var_intersect_j.object;

    if ( tmp_binop_right_20 == NULL )
    {
        Py_DECREF( tmp_called_9 );
        Py_DECREF( tmp_call_arg_element_21 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 963 ], 57, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 88;
        goto try_finally_handler_1;
    }

    tmp_binop_left_19 = BINARY_OPERATION_ADD( tmp_binop_left_20, tmp_binop_right_20 );
    if ( tmp_binop_left_19 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_called_9 );
        Py_DECREF( tmp_call_arg_element_21 );

        frame_function->f_lineno = 88;
        goto try_finally_handler_1;
    }
    tmp_binop_right_19 = var_i.object;

    if ( tmp_binop_right_19 == NULL )
    {
        Py_DECREF( tmp_called_9 );
        Py_DECREF( tmp_call_arg_element_21 );
        Py_DECREF( tmp_binop_left_19 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 182 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 89;
        goto try_finally_handler_1;
    }

    tmp_list_element_2 = BINARY_OPERATION_ADD( tmp_binop_left_19, tmp_binop_right_19 );
    Py_DECREF( tmp_binop_left_19 );
    if ( tmp_list_element_2 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_called_9 );
        Py_DECREF( tmp_call_arg_element_21 );

        frame_function->f_lineno = 88;
        goto try_finally_handler_1;
    }
    PyList_SET_ITEM( tmp_call_arg_element_21, 0, tmp_list_element_2 );
    tmp_binop_left_21 = par_current_index.object;

    if ( tmp_binop_left_21 == NULL )
    {
        Py_DECREF( tmp_called_9 );
        Py_DECREF( tmp_call_arg_element_21 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1221 ], 59, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 90;
        goto try_finally_handler_1;
    }

    tmp_binop_right_21 = var_intersect_i.object;

    if ( tmp_binop_right_21 == NULL )
    {
        Py_DECREF( tmp_called_9 );
        Py_DECREF( tmp_call_arg_element_21 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 857 ], 57, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 90;
        goto try_finally_handler_1;
    }

    tmp_list_element_2 = BINARY_OPERATION_ADD( tmp_binop_left_21, tmp_binop_right_21 );
    if ( tmp_list_element_2 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_called_9 );
        Py_DECREF( tmp_call_arg_element_21 );

        frame_function->f_lineno = 90;
        goto try_finally_handler_1;
    }
    PyList_SET_ITEM( tmp_call_arg_element_21, 1, tmp_list_element_2 );
    tmp_binop_left_22 = const_float__minus_1_0;
    tmp_binop_right_22 = var_crossing_sign.object;

    tmp_list_element_2 = BINARY_OPERATION_MUL( tmp_binop_left_22, tmp_binop_right_22 );
    if ( tmp_list_element_2 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_called_9 );
        Py_DECREF( tmp_call_arg_element_21 );

        frame_function->f_lineno = 91;
        goto try_finally_handler_1;
    }
    PyList_SET_ITEM( tmp_call_arg_element_21, 2, tmp_list_element_2 );
    tmp_binop_left_23 = var_crossing_sign.object;

    tmp_binop_right_23 = var_crossing_direction.object;

    tmp_list_element_2 = BINARY_OPERATION_MUL( tmp_binop_left_23, tmp_binop_right_23 );
    if ( tmp_list_element_2 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_called_9 );
        Py_DECREF( tmp_call_arg_element_21 );

        frame_function->f_lineno = 92;
        goto try_finally_handler_1;
    }
    PyList_SET_ITEM( tmp_call_arg_element_21, 3, tmp_list_element_2 );
    frame_function->f_lineno = 92;
    tmp_unused = CALL_FUNCTION_WITH_ARGS1( tmp_called_9, tmp_call_arg_element_21 );
    Py_DECREF( tmp_called_9 );
    Py_DECREF( tmp_call_arg_element_21 );
    if ( tmp_unused == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 92;
        goto try_finally_handler_1;
    }
    Py_DECREF( tmp_unused );
    branch_no_3:;
    tmp_assign_source_30 = var_i.object;

    if ( tmp_assign_source_30 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 182 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 93;
        goto try_finally_handler_1;
    }

    if (tmp_inplace_assign_1__inplace_start.object == NULL)
    {
        tmp_inplace_assign_1__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_30 );
    }
    else
    {
        PyObject *old = tmp_inplace_assign_1__inplace_start.object;
        tmp_inplace_assign_1__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_30 );
        Py_DECREF( old );
    }
    // Tried code
    tmp_binop_left_24 = tmp_inplace_assign_1__inplace_start.object;

    tmp_binop_right_24 = const_int_pos_1;
    tmp_assign_source_31 = BINARY_OPERATION( PyNumber_InPlaceAdd, tmp_binop_left_24, tmp_binop_right_24 );
    if ( tmp_assign_source_31 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 93;
        goto try_finally_handler_5;
    }
    if (tmp_inplace_assign_1__inplace_end.object == NULL)
    {
        tmp_inplace_assign_1__inplace_end.object = tmp_assign_source_31;
    }
    else
    {
        PyObject *old = tmp_inplace_assign_1__inplace_end.object;
        tmp_inplace_assign_1__inplace_end.object = tmp_assign_source_31;
        Py_DECREF( old );
    }
    tmp_compare_left_2 = tmp_inplace_assign_1__inplace_start.object;

    tmp_compare_right_2 = tmp_inplace_assign_1__inplace_end.object;

    tmp_isnot_1 = ( tmp_compare_left_2 != tmp_compare_right_2 );
    if (tmp_isnot_1)
    {
        goto branch_yes_4;
    }
    else
    {
        goto branch_no_4;
    }
    branch_yes_4:;
    tmp_assign_source_32 = tmp_inplace_assign_1__inplace_end.object;

    if (var_i.object == NULL)
    {
        var_i.object = INCREASE_REFCOUNT( tmp_assign_source_32 );
    }
    else
    {
        PyObject *old = var_i.object;
        var_i.object = INCREASE_REFCOUNT( tmp_assign_source_32 );
        Py_DECREF( old );
    }
    branch_no_4:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_5:;
    exception_keeper_type_4 = exception_type;
    exception_keeper_value_4 = exception_value;
    exception_keeper_tb_4 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    tmp_tried_lineno_2 = frame_function->f_lineno;
    tmp_result = tmp_inplace_assign_1__inplace_start.object != NULL;
    if ( tmp_result == true )
    {
        Py_DECREF( tmp_inplace_assign_1__inplace_start.object );
        tmp_inplace_assign_1__inplace_start.object = NULL;
    }

    assert( tmp_result != false );
    Py_XDECREF( tmp_inplace_assign_1__inplace_end.object );
    tmp_inplace_assign_1__inplace_end.object = NULL;

    frame_function->f_lineno = tmp_tried_lineno_2;
    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_4 != NULL )
    {
        exception_type = exception_keeper_type_4;
        exception_value = exception_keeper_value_4;
        exception_tb = exception_keeper_tb_4;

        goto try_finally_handler_1;
    }

    goto finally_end_4;
    finally_end_4:;
    goto branch_end_2;
    branch_no_2:;
    tmp_compare_left_3 = par_jump_mode.object;

    if ( tmp_compare_left_3 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1342 ], 55, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 95;
        goto try_finally_handler_1;
    }

    tmp_compare_right_3 = const_int_pos_3;
    tmp_cmp_Eq_1 = RICH_COMPARE_BOOL_EQ( tmp_compare_left_3, tmp_compare_right_3 );
    if ( tmp_cmp_Eq_1 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 95;
        goto try_finally_handler_1;
    }
    if (tmp_cmp_Eq_1 == 1)
    {
        goto branch_yes_5;
    }
    else
    {
        goto branch_no_5;
    }
    branch_yes_5:;
    tmp_assign_source_33 = var_i.object;

    if ( tmp_assign_source_33 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 182 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 96;
        goto try_finally_handler_1;
    }

    if (tmp_inplace_assign_2__inplace_start.object == NULL)
    {
        tmp_inplace_assign_2__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_33 );
    }
    else
    {
        PyObject *old = tmp_inplace_assign_2__inplace_start.object;
        tmp_inplace_assign_2__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_33 );
        Py_DECREF( old );
    }
    // Tried code
    tmp_binop_left_25 = tmp_inplace_assign_2__inplace_start.object;

    tmp_binop_right_25 = const_int_pos_1;
    tmp_assign_source_34 = BINARY_OPERATION( PyNumber_InPlaceAdd, tmp_binop_left_25, tmp_binop_right_25 );
    if ( tmp_assign_source_34 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 96;
        goto try_finally_handler_6;
    }
    if (tmp_inplace_assign_2__inplace_end.object == NULL)
    {
        tmp_inplace_assign_2__inplace_end.object = tmp_assign_source_34;
    }
    else
    {
        PyObject *old = tmp_inplace_assign_2__inplace_end.object;
        tmp_inplace_assign_2__inplace_end.object = tmp_assign_source_34;
        Py_DECREF( old );
    }
    tmp_compare_left_4 = tmp_inplace_assign_2__inplace_start.object;

    tmp_compare_right_4 = tmp_inplace_assign_2__inplace_end.object;

    tmp_isnot_2 = ( tmp_compare_left_4 != tmp_compare_right_4 );
    if (tmp_isnot_2)
    {
        goto branch_yes_6;
    }
    else
    {
        goto branch_no_6;
    }
    branch_yes_6:;
    tmp_assign_source_35 = tmp_inplace_assign_2__inplace_end.object;

    if (var_i.object == NULL)
    {
        var_i.object = INCREASE_REFCOUNT( tmp_assign_source_35 );
    }
    else
    {
        PyObject *old = var_i.object;
        var_i.object = INCREASE_REFCOUNT( tmp_assign_source_35 );
        Py_DECREF( old );
    }
    branch_no_6:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_6:;
    exception_keeper_type_5 = exception_type;
    exception_keeper_value_5 = exception_value;
    exception_keeper_tb_5 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    tmp_tried_lineno_3 = frame_function->f_lineno;
    tmp_result = tmp_inplace_assign_2__inplace_start.object != NULL;
    if ( tmp_result == true )
    {
        Py_DECREF( tmp_inplace_assign_2__inplace_start.object );
        tmp_inplace_assign_2__inplace_start.object = NULL;
    }

    assert( tmp_result != false );
    Py_XDECREF( tmp_inplace_assign_2__inplace_end.object );
    tmp_inplace_assign_2__inplace_end.object = NULL;

    frame_function->f_lineno = tmp_tried_lineno_3;
    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_5 != NULL )
    {
        exception_type = exception_keeper_type_5;
        exception_value = exception_keeper_value_5;
        exception_tb = exception_keeper_tb_5;

        goto try_finally_handler_1;
    }

    goto finally_end_5;
    finally_end_5:;
    tmp_assign_source_36 = const_int_pos_1;
    if (var_already_jumped.object == NULL)
    {
        var_already_jumped.object = INCREASE_REFCOUNT( tmp_assign_source_36 );
    }
    else
    {
        PyObject *old = var_already_jumped.object;
        var_already_jumped.object = INCREASE_REFCOUNT( tmp_assign_source_36 );
        Py_DECREF( old );
    }
    goto branch_end_5;
    branch_no_5:;
    tmp_compare_left_5 = par_jump_mode.object;

    if ( tmp_compare_left_5 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1342 ], 55, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 98;
        goto try_finally_handler_1;
    }

    tmp_compare_right_5 = const_int_pos_2;
    tmp_cmp_Eq_2 = RICH_COMPARE_BOOL_EQ( tmp_compare_left_5, tmp_compare_right_5 );
    if ( tmp_cmp_Eq_2 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 98;
        goto try_finally_handler_1;
    }
    if (tmp_cmp_Eq_2 == 1)
    {
        goto branch_yes_7;
    }
    else
    {
        goto branch_no_7;
    }
    branch_yes_7:;
    tmp_called_10 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_floor );

    if (unlikely( tmp_called_10 == NULL ))
    {
        tmp_called_10 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_floor );
    }

    if ( tmp_called_10 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1397 ], 34, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 99;
        goto try_finally_handler_1;
    }

    tmp_binop_left_27 = var_distance.object;

    if ( tmp_binop_left_27 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1431 ], 54, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 99;
        goto try_finally_handler_1;
    }

    tmp_binop_right_27 = par_max_segment_length.object;

    if ( tmp_binop_right_27 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 118 ], 64, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 99;
        goto try_finally_handler_1;
    }

    tmp_call_arg_element_22 = BINARY_OPERATION_DIV( tmp_binop_left_27, tmp_binop_right_27 );
    if ( tmp_call_arg_element_22 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 99;
        goto try_finally_handler_1;
    }
    frame_function->f_lineno = 99;
    tmp_binop_left_26 = CALL_FUNCTION_WITH_ARGS1( tmp_called_10, tmp_call_arg_element_22 );
    Py_DECREF( tmp_call_arg_element_22 );
    if ( tmp_binop_left_26 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 99;
        goto try_finally_handler_1;
    }
    tmp_binop_right_26 = const_int_pos_1;
    tmp_assign_source_37 = BINARY_OPERATION_SUB( tmp_binop_left_26, tmp_binop_right_26 );
    Py_DECREF( tmp_binop_left_26 );
    if ( tmp_assign_source_37 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 99;
        goto try_finally_handler_1;
    }
    if (var_num_jumps.object == NULL)
    {
        var_num_jumps.object = tmp_assign_source_37;
    }
    else
    {
        PyObject *old = var_num_jumps.object;
        var_num_jumps.object = tmp_assign_source_37;
        Py_DECREF( old );
    }
    tmp_compare_left_6 = var_num_jumps.object;

    tmp_compare_right_6 = const_int_pos_1;
    tmp_cmp_Lt_2 = RICH_COMPARE_BOOL_LT( tmp_compare_left_6, tmp_compare_right_6 );
    if ( tmp_cmp_Lt_2 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 100;
        goto try_finally_handler_1;
    }
    if (tmp_cmp_Lt_2 == 1)
    {
        goto branch_yes_8;
    }
    else
    {
        goto branch_no_8;
    }
    branch_yes_8:;
    tmp_assign_source_38 = const_int_pos_1;
    assert( var_num_jumps.object != NULL );
    {
        PyObject *old = var_num_jumps.object;
        var_num_jumps.object = INCREASE_REFCOUNT( tmp_assign_source_38 );
        Py_DECREF( old );
    }

    branch_no_8:;
    tmp_assign_source_39 = var_i.object;

    if ( tmp_assign_source_39 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 182 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 102;
        goto try_finally_handler_1;
    }

    if (tmp_inplace_assign_3__inplace_start.object == NULL)
    {
        tmp_inplace_assign_3__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_39 );
    }
    else
    {
        PyObject *old = tmp_inplace_assign_3__inplace_start.object;
        tmp_inplace_assign_3__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_39 );
        Py_DECREF( old );
    }
    // Tried code
    tmp_binop_left_28 = tmp_inplace_assign_3__inplace_start.object;

    tmp_binop_right_28 = var_num_jumps.object;

    tmp_assign_source_40 = BINARY_OPERATION( PyNumber_InPlaceAdd, tmp_binop_left_28, tmp_binop_right_28 );
    if ( tmp_assign_source_40 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 102;
        goto try_finally_handler_7;
    }
    if (tmp_inplace_assign_3__inplace_end.object == NULL)
    {
        tmp_inplace_assign_3__inplace_end.object = tmp_assign_source_40;
    }
    else
    {
        PyObject *old = tmp_inplace_assign_3__inplace_end.object;
        tmp_inplace_assign_3__inplace_end.object = tmp_assign_source_40;
        Py_DECREF( old );
    }
    tmp_compare_left_7 = tmp_inplace_assign_3__inplace_start.object;

    tmp_compare_right_7 = tmp_inplace_assign_3__inplace_end.object;

    tmp_isnot_3 = ( tmp_compare_left_7 != tmp_compare_right_7 );
    if (tmp_isnot_3)
    {
        goto branch_yes_9;
    }
    else
    {
        goto branch_no_9;
    }
    branch_yes_9:;
    tmp_assign_source_41 = tmp_inplace_assign_3__inplace_end.object;

    if (var_i.object == NULL)
    {
        var_i.object = INCREASE_REFCOUNT( tmp_assign_source_41 );
    }
    else
    {
        PyObject *old = var_i.object;
        var_i.object = INCREASE_REFCOUNT( tmp_assign_source_41 );
        Py_DECREF( old );
    }
    branch_no_9:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_7:;
    exception_keeper_type_6 = exception_type;
    exception_keeper_value_6 = exception_value;
    exception_keeper_tb_6 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    tmp_tried_lineno_4 = frame_function->f_lineno;
    tmp_result = tmp_inplace_assign_3__inplace_start.object != NULL;
    if ( tmp_result == true )
    {
        Py_DECREF( tmp_inplace_assign_3__inplace_start.object );
        tmp_inplace_assign_3__inplace_start.object = NULL;
    }

    assert( tmp_result != false );
    Py_XDECREF( tmp_inplace_assign_3__inplace_end.object );
    tmp_inplace_assign_3__inplace_end.object = NULL;

    frame_function->f_lineno = tmp_tried_lineno_4;
    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_6 != NULL )
    {
        exception_type = exception_keeper_type_6;
        exception_value = exception_keeper_value_6;
        exception_tb = exception_keeper_tb_6;

        goto try_finally_handler_1;
    }

    goto finally_end_6;
    finally_end_6:;
    tmp_assign_source_42 = const_int_pos_1;
    if (var_already_jumped.object == NULL)
    {
        var_already_jumped.object = INCREASE_REFCOUNT( tmp_assign_source_42 );
    }
    else
    {
        PyObject *old = var_already_jumped.object;
        var_already_jumped.object = INCREASE_REFCOUNT( tmp_assign_source_42 );
        Py_DECREF( old );
    }
    goto branch_end_7;
    branch_no_7:;
    tmp_assign_source_43 = const_float_0_0;
    if (var_distance_travelled.object == NULL)
    {
        var_distance_travelled.object = INCREASE_REFCOUNT( tmp_assign_source_43 );
    }
    else
    {
        PyObject *old = var_distance_travelled.object;
        var_distance_travelled.object = INCREASE_REFCOUNT( tmp_assign_source_43 );
        Py_DECREF( old );
    }
    tmp_assign_source_44 = const_int_0;
    if (var_jumps.object == NULL)
    {
        var_jumps.object = INCREASE_REFCOUNT( tmp_assign_source_44 );
    }
    else
    {
        PyObject *old = var_jumps.object;
        var_jumps.object = INCREASE_REFCOUNT( tmp_assign_source_44 );
        Py_DECREF( old );
    }
    loop_start_2:;
    tmp_break_1 = false;
    // Tried code
    tmp_cond_value_4 = NULL;
    // Tried code
    tmp_compexpr_left_2 = var_distance_travelled.object;

    if ( tmp_compexpr_left_2 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1485 ], 64, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 107;
        goto try_finally_handler_9;
    }

    tmp_binop_left_29 = var_distance.object;

    if ( tmp_binop_left_29 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1431 ], 54, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 107;
        goto try_finally_handler_9;
    }

    tmp_binop_right_29 = par_max_segment_length.object;

    if ( tmp_binop_right_29 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 118 ], 64, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 107;
        goto try_finally_handler_9;
    }

    tmp_compexpr_right_2 = BINARY_OPERATION_SUB( tmp_binop_left_29, tmp_binop_right_29 );
    if ( tmp_compexpr_right_2 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 107;
        goto try_finally_handler_9;
    }
    tmp_assign_source_45 = RICH_COMPARE_LT( tmp_compexpr_left_2, tmp_compexpr_right_2 );
    Py_DECREF( tmp_compexpr_right_2 );
    if ( tmp_assign_source_45 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 107;
        goto try_finally_handler_9;
    }
    if (tmp_and_1__value_1.object == NULL)
    {
        tmp_and_1__value_1.object = tmp_assign_source_45;
    }
    else
    {
        PyObject *old = tmp_and_1__value_1.object;
        tmp_and_1__value_1.object = tmp_assign_source_45;
        Py_DECREF( old );
    }
    tmp_cond_value_5 = tmp_and_1__value_1.object;

    tmp_cond_truth_5 = CHECK_IF_TRUE( tmp_cond_value_5 );
    if ( tmp_cond_truth_5 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 107;
        goto try_finally_handler_9;
    }
    if (tmp_cond_truth_5 == 1)
    {
        goto condexpr_true_2;
    }
    else
    {
        goto condexpr_false_2;
    }
    condexpr_true_2:;
    tmp_cond_value_4 = NULL;
    // Tried code
    tmp_result = tmp_and_1__value_1.object != NULL;
    if ( tmp_result == true )
    {
        Py_DECREF( tmp_and_1__value_1.object );
        tmp_and_1__value_1.object = NULL;
    }

    assert( tmp_result != false );
    tmp_compexpr_left_3 = var_i.object;

    if ( tmp_compexpr_left_3 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 182 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 108;
        goto try_finally_handler_10;
    }

    tmp_len_arg_2 = par_points.object;

    if ( tmp_len_arg_2 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 229 ], 52, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 108;
        goto try_finally_handler_10;
    }

    tmp_compexpr_right_3 = BUILTIN_LEN( tmp_len_arg_2 );
    if ( tmp_compexpr_right_3 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 108;
        goto try_finally_handler_10;
    }
    tmp_cond_value_4 = RICH_COMPARE_LT( tmp_compexpr_left_3, tmp_compexpr_right_3 );
    Py_DECREF( tmp_compexpr_right_3 );
    if ( tmp_cond_value_4 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 108;
        goto try_finally_handler_10;
    }
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_10:;
    exception_keeper_type_7 = exception_type;
    exception_keeper_value_7 = exception_value;
    exception_keeper_tb_7 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_7 != NULL )
    {
        exception_type = exception_keeper_type_7;
        exception_value = exception_keeper_value_7;
        exception_tb = exception_keeper_tb_7;

        goto try_finally_handler_9;
    }

    goto finally_end_7;
    finally_end_7:;
    goto condexpr_end_2;
    condexpr_false_2:;
    tmp_cond_value_4 = tmp_and_1__value_1.object;

    Py_INCREF( tmp_cond_value_4 );
    condexpr_end_2:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_9:;
    exception_keeper_type_8 = exception_type;
    exception_keeper_value_8 = exception_value;
    exception_keeper_tb_8 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_8 != NULL )
    {
        exception_type = exception_keeper_type_8;
        exception_value = exception_keeper_value_8;
        exception_tb = exception_keeper_tb_8;

        goto try_finally_handler_8;
    }

    goto finally_end_8;
    finally_end_8:;
    tmp_cond_truth_4 = CHECK_IF_TRUE( tmp_cond_value_4 );
    if ( tmp_cond_truth_4 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_cond_value_4 );

        frame_function->f_lineno = 107;
        goto try_finally_handler_8;
    }
    Py_DECREF( tmp_cond_value_4 );
    if (tmp_cond_truth_4 == 1)
    {
        goto branch_no_10;
    }
    else
    {
        goto branch_yes_10;
    }
    branch_yes_10:;
    tmp_break_1 = true;
    goto try_finally_handler_start_1;
    branch_no_10:;
    try_finally_handler_start_1:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_8:;
    exception_keeper_type_9 = exception_type;
    exception_keeper_value_9 = exception_value;
    exception_keeper_tb_9 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    tmp_tried_lineno_5 = frame_function->f_lineno;
    Py_XDECREF( tmp_and_1__value_1.object );
    tmp_and_1__value_1.object = NULL;

    frame_function->f_lineno = tmp_tried_lineno_5;
    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_9 != NULL )
    {
        exception_type = exception_keeper_type_9;
        exception_value = exception_keeper_value_9;
        exception_tb = exception_keeper_tb_9;

        goto try_finally_handler_1;
    }

    // Break if entered via break.
    if ( tmp_break_1 )
    {

    goto loop_end_2;
    }
    goto finally_end_9;
    finally_end_9:;
    tmp_assign_source_46 = var_jumps.object;

    if ( tmp_assign_source_46 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1549 ], 51, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 109;
        goto try_finally_handler_1;
    }

    if (tmp_inplace_assign_4__inplace_start.object == NULL)
    {
        tmp_inplace_assign_4__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_46 );
    }
    else
    {
        PyObject *old = tmp_inplace_assign_4__inplace_start.object;
        tmp_inplace_assign_4__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_46 );
        Py_DECREF( old );
    }
    // Tried code
    tmp_binop_left_30 = tmp_inplace_assign_4__inplace_start.object;

    tmp_binop_right_30 = const_int_pos_1;
    tmp_assign_source_47 = BINARY_OPERATION( PyNumber_InPlaceAdd, tmp_binop_left_30, tmp_binop_right_30 );
    if ( tmp_assign_source_47 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 109;
        goto try_finally_handler_11;
    }
    if (tmp_inplace_assign_4__inplace_end.object == NULL)
    {
        tmp_inplace_assign_4__inplace_end.object = tmp_assign_source_47;
    }
    else
    {
        PyObject *old = tmp_inplace_assign_4__inplace_end.object;
        tmp_inplace_assign_4__inplace_end.object = tmp_assign_source_47;
        Py_DECREF( old );
    }
    tmp_compare_left_8 = tmp_inplace_assign_4__inplace_start.object;

    tmp_compare_right_8 = tmp_inplace_assign_4__inplace_end.object;

    tmp_isnot_4 = ( tmp_compare_left_8 != tmp_compare_right_8 );
    if (tmp_isnot_4)
    {
        goto branch_yes_11;
    }
    else
    {
        goto branch_no_11;
    }
    branch_yes_11:;
    tmp_assign_source_48 = tmp_inplace_assign_4__inplace_end.object;

    if (var_jumps.object == NULL)
    {
        var_jumps.object = INCREASE_REFCOUNT( tmp_assign_source_48 );
    }
    else
    {
        PyObject *old = var_jumps.object;
        var_jumps.object = INCREASE_REFCOUNT( tmp_assign_source_48 );
        Py_DECREF( old );
    }
    branch_no_11:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_11:;
    exception_keeper_type_10 = exception_type;
    exception_keeper_value_10 = exception_value;
    exception_keeper_tb_10 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    tmp_tried_lineno_6 = frame_function->f_lineno;
    tmp_result = tmp_inplace_assign_4__inplace_start.object != NULL;
    if ( tmp_result == true )
    {
        Py_DECREF( tmp_inplace_assign_4__inplace_start.object );
        tmp_inplace_assign_4__inplace_start.object = NULL;
    }

    assert( tmp_result != false );
    Py_XDECREF( tmp_inplace_assign_4__inplace_end.object );
    tmp_inplace_assign_4__inplace_end.object = NULL;

    frame_function->f_lineno = tmp_tried_lineno_6;
    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_10 != NULL )
    {
        exception_type = exception_keeper_type_10;
        exception_value = exception_keeper_value_10;
        exception_tb = exception_keeper_tb_10;

        goto try_finally_handler_1;
    }

    goto finally_end_10;
    finally_end_10:;
    tmp_assign_source_49 = var_distance_travelled.object;

    if ( tmp_assign_source_49 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1485 ], 64, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 110;
        goto try_finally_handler_1;
    }

    if (tmp_inplace_assign_5__inplace_start.object == NULL)
    {
        tmp_inplace_assign_5__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_49 );
    }
    else
    {
        PyObject *old = tmp_inplace_assign_5__inplace_start.object;
        tmp_inplace_assign_5__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_49 );
        Py_DECREF( old );
    }
    // Tried code
    tmp_binop_left_31 = tmp_inplace_assign_5__inplace_start.object;

    tmp_subscr_target_20 = par_segment_lengths.object;

    if ( tmp_subscr_target_20 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1600 ], 61, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 110;
        goto try_finally_handler_12;
    }

    tmp_subscr_subscript_20 = var_i.object;

    if ( tmp_subscr_subscript_20 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 182 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 110;
        goto try_finally_handler_12;
    }

    tmp_binop_right_31 = LOOKUP_SUBSCRIPT( tmp_subscr_target_20, tmp_subscr_subscript_20 );
    if ( tmp_binop_right_31 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 110;
        goto try_finally_handler_12;
    }
    tmp_assign_source_50 = BINARY_OPERATION( PyNumber_InPlaceAdd, tmp_binop_left_31, tmp_binop_right_31 );
    Py_DECREF( tmp_binop_right_31 );
    if ( tmp_assign_source_50 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 110;
        goto try_finally_handler_12;
    }
    if (tmp_inplace_assign_5__inplace_end.object == NULL)
    {
        tmp_inplace_assign_5__inplace_end.object = tmp_assign_source_50;
    }
    else
    {
        PyObject *old = tmp_inplace_assign_5__inplace_end.object;
        tmp_inplace_assign_5__inplace_end.object = tmp_assign_source_50;
        Py_DECREF( old );
    }
    tmp_compare_left_9 = tmp_inplace_assign_5__inplace_start.object;

    tmp_compare_right_9 = tmp_inplace_assign_5__inplace_end.object;

    tmp_isnot_5 = ( tmp_compare_left_9 != tmp_compare_right_9 );
    if (tmp_isnot_5)
    {
        goto branch_yes_12;
    }
    else
    {
        goto branch_no_12;
    }
    branch_yes_12:;
    tmp_assign_source_51 = tmp_inplace_assign_5__inplace_end.object;

    if (var_distance_travelled.object == NULL)
    {
        var_distance_travelled.object = INCREASE_REFCOUNT( tmp_assign_source_51 );
    }
    else
    {
        PyObject *old = var_distance_travelled.object;
        var_distance_travelled.object = INCREASE_REFCOUNT( tmp_assign_source_51 );
        Py_DECREF( old );
    }
    branch_no_12:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_12:;
    exception_keeper_type_11 = exception_type;
    exception_keeper_value_11 = exception_value;
    exception_keeper_tb_11 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    tmp_tried_lineno_7 = frame_function->f_lineno;
    tmp_result = tmp_inplace_assign_5__inplace_start.object != NULL;
    if ( tmp_result == true )
    {
        Py_DECREF( tmp_inplace_assign_5__inplace_start.object );
        tmp_inplace_assign_5__inplace_start.object = NULL;
    }

    assert( tmp_result != false );
    Py_XDECREF( tmp_inplace_assign_5__inplace_end.object );
    tmp_inplace_assign_5__inplace_end.object = NULL;

    frame_function->f_lineno = tmp_tried_lineno_7;
    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_11 != NULL )
    {
        exception_type = exception_keeper_type_11;
        exception_value = exception_keeper_value_11;
        exception_tb = exception_keeper_tb_11;

        goto try_finally_handler_1;
    }

    goto finally_end_11;
    finally_end_11:;
    tmp_assign_source_52 = var_i.object;

    if ( tmp_assign_source_52 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 182 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 111;
        goto try_finally_handler_1;
    }

    if (tmp_inplace_assign_6__inplace_start.object == NULL)
    {
        tmp_inplace_assign_6__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_52 );
    }
    else
    {
        PyObject *old = tmp_inplace_assign_6__inplace_start.object;
        tmp_inplace_assign_6__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_52 );
        Py_DECREF( old );
    }
    // Tried code
    tmp_binop_left_32 = tmp_inplace_assign_6__inplace_start.object;

    tmp_binop_right_32 = const_int_pos_1;
    tmp_assign_source_53 = BINARY_OPERATION( PyNumber_InPlaceAdd, tmp_binop_left_32, tmp_binop_right_32 );
    if ( tmp_assign_source_53 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 111;
        goto try_finally_handler_13;
    }
    if (tmp_inplace_assign_6__inplace_end.object == NULL)
    {
        tmp_inplace_assign_6__inplace_end.object = tmp_assign_source_53;
    }
    else
    {
        PyObject *old = tmp_inplace_assign_6__inplace_end.object;
        tmp_inplace_assign_6__inplace_end.object = tmp_assign_source_53;
        Py_DECREF( old );
    }
    tmp_compare_left_10 = tmp_inplace_assign_6__inplace_start.object;

    tmp_compare_right_10 = tmp_inplace_assign_6__inplace_end.object;

    tmp_isnot_6 = ( tmp_compare_left_10 != tmp_compare_right_10 );
    if (tmp_isnot_6)
    {
        goto branch_yes_13;
    }
    else
    {
        goto branch_no_13;
    }
    branch_yes_13:;
    tmp_assign_source_54 = tmp_inplace_assign_6__inplace_end.object;

    if (var_i.object == NULL)
    {
        var_i.object = INCREASE_REFCOUNT( tmp_assign_source_54 );
    }
    else
    {
        PyObject *old = var_i.object;
        var_i.object = INCREASE_REFCOUNT( tmp_assign_source_54 );
        Py_DECREF( old );
    }
    branch_no_13:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_13:;
    exception_keeper_type_12 = exception_type;
    exception_keeper_value_12 = exception_value;
    exception_keeper_tb_12 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    tmp_tried_lineno_8 = frame_function->f_lineno;
    tmp_result = tmp_inplace_assign_6__inplace_start.object != NULL;
    if ( tmp_result == true )
    {
        Py_DECREF( tmp_inplace_assign_6__inplace_start.object );
        tmp_inplace_assign_6__inplace_start.object = NULL;
    }

    assert( tmp_result != false );
    Py_XDECREF( tmp_inplace_assign_6__inplace_end.object );
    tmp_inplace_assign_6__inplace_end.object = NULL;

    frame_function->f_lineno = tmp_tried_lineno_8;
    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_12 != NULL )
    {
        exception_type = exception_keeper_type_12;
        exception_value = exception_keeper_value_12;
        exception_tb = exception_keeper_tb_12;

        goto try_finally_handler_1;
    }

    goto finally_end_12;
    finally_end_12:;
    if ( CONSIDER_THREADING() == false )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 107;
        goto try_finally_handler_1;
    }
    goto loop_start_2;
    loop_end_2:;
    tmp_compare_left_11 = var_jumps.object;

    if ( tmp_compare_left_11 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1549 ], 51, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 112;
        goto try_finally_handler_1;
    }

    tmp_compare_right_11 = const_int_pos_1;
    tmp_cmp_Gt_1 = RICH_COMPARE_BOOL_GT( tmp_compare_left_11, tmp_compare_right_11 );
    if ( tmp_cmp_Gt_1 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 112;
        goto try_finally_handler_1;
    }
    if (tmp_cmp_Gt_1 == 1)
    {
        goto branch_yes_14;
    }
    else
    {
        goto branch_no_14;
    }
    branch_yes_14:;
    tmp_assign_source_55 = var_i.object;

    if ( tmp_assign_source_55 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 182 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 113;
        goto try_finally_handler_1;
    }

    if (tmp_inplace_assign_7__inplace_start.object == NULL)
    {
        tmp_inplace_assign_7__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_55 );
    }
    else
    {
        PyObject *old = tmp_inplace_assign_7__inplace_start.object;
        tmp_inplace_assign_7__inplace_start.object = INCREASE_REFCOUNT( tmp_assign_source_55 );
        Py_DECREF( old );
    }
    // Tried code
    tmp_binop_left_33 = tmp_inplace_assign_7__inplace_start.object;

    tmp_binop_right_33 = const_int_pos_2;
    tmp_assign_source_56 = BINARY_OPERATION( PyNumber_InPlaceSubtract, tmp_binop_left_33, tmp_binop_right_33 );
    if ( tmp_assign_source_56 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 113;
        goto try_finally_handler_14;
    }
    if (tmp_inplace_assign_7__inplace_end.object == NULL)
    {
        tmp_inplace_assign_7__inplace_end.object = tmp_assign_source_56;
    }
    else
    {
        PyObject *old = tmp_inplace_assign_7__inplace_end.object;
        tmp_inplace_assign_7__inplace_end.object = tmp_assign_source_56;
        Py_DECREF( old );
    }
    tmp_compare_left_12 = tmp_inplace_assign_7__inplace_start.object;

    tmp_compare_right_12 = tmp_inplace_assign_7__inplace_end.object;

    tmp_isnot_7 = ( tmp_compare_left_12 != tmp_compare_right_12 );
    if (tmp_isnot_7)
    {
        goto branch_yes_15;
    }
    else
    {
        goto branch_no_15;
    }
    branch_yes_15:;
    tmp_assign_source_57 = tmp_inplace_assign_7__inplace_end.object;

    if (var_i.object == NULL)
    {
        var_i.object = INCREASE_REFCOUNT( tmp_assign_source_57 );
    }
    else
    {
        PyObject *old = var_i.object;
        var_i.object = INCREASE_REFCOUNT( tmp_assign_source_57 );
        Py_DECREF( old );
    }
    branch_no_15:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_14:;
    exception_keeper_type_13 = exception_type;
    exception_keeper_value_13 = exception_value;
    exception_keeper_tb_13 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    tmp_tried_lineno_9 = frame_function->f_lineno;
    tmp_result = tmp_inplace_assign_7__inplace_start.object != NULL;
    if ( tmp_result == true )
    {
        Py_DECREF( tmp_inplace_assign_7__inplace_start.object );
        tmp_inplace_assign_7__inplace_start.object = NULL;
    }

    assert( tmp_result != false );
    Py_XDECREF( tmp_inplace_assign_7__inplace_end.object );
    tmp_inplace_assign_7__inplace_end.object = NULL;

    frame_function->f_lineno = tmp_tried_lineno_9;
    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_13 != NULL )
    {
        exception_type = exception_keeper_type_13;
        exception_value = exception_keeper_value_13;
        exception_tb = exception_keeper_tb_13;

        goto try_finally_handler_1;
    }

    goto finally_end_13;
    finally_end_13:;
    branch_no_14:;
    tmp_assign_source_58 = const_int_pos_1;
    if (var_already_jumped.object == NULL)
    {
        var_already_jumped.object = INCREASE_REFCOUNT( tmp_assign_source_58 );
    }
    else
    {
        PyObject *old = var_already_jumped.object;
        var_already_jumped.object = INCREASE_REFCOUNT( tmp_assign_source_58 );
        Py_DECREF( old );
    }
    branch_end_7:;
    branch_end_5:;
    branch_end_2:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_1:;
    exception_keeper_type_14 = exception_type;
    exception_keeper_value_14 = exception_value;
    exception_keeper_tb_14 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    tmp_tried_lineno_10 = frame_function->f_lineno;
    Py_XDECREF( tmp_or_1__value_1.object );
    tmp_or_1__value_1.object = NULL;

    frame_function->f_lineno = tmp_tried_lineno_10;
    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_14 != NULL )
    {
        exception_type = exception_keeper_type_14;
        exception_value = exception_keeper_value_14;
        exception_tb = exception_keeper_tb_14;

        goto frame_exception_exit_1;
    }

    goto finally_end_14;
    finally_end_14:;
    if ( CONSIDER_THREADING() == false )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 58;
        goto frame_exception_exit_1;
    }
    goto loop_start_1;
    loop_end_1:;
    tmp_return_value = var_crossings.object;

    if ( tmp_return_value == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1166 ], 55, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 118;
        goto frame_exception_exit_1;
    }

    Py_INCREF( tmp_return_value );
    goto frame_return_exit_1;

#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif
    // Put the previous frame back on top.
    popFrameStack();
#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );
    goto frame_no_exception_1;
    frame_return_exit_1:;
#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif
    popFrameStack();
#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );
    goto function_return_exit;
    frame_exception_exit_1:;
#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif

    if ( exception_tb == NULL )
    {
        exception_tb = MAKE_TRACEBACK( INCREASE_REFCOUNT( frame_function ) );
    }
    else if ( exception_tb->tb_frame != frame_function )
    {
        PyTracebackObject *traceback_new = (PyTracebackObject *)MAKE_TRACEBACK( INCREASE_REFCOUNT( frame_function ) );
        traceback_new->tb_next = exception_tb;
        exception_tb = traceback_new;
    }


    tmp_frame_locals = PyDict_New();
    if ((var_crossings.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_crossings,
            var_crossings.object
        );

    }
    if ((var_vx.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_vx,
            var_vx.object
        );

    }
    if ((var_vy.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_vy,
            var_vy.object
        );

    }
    if ((var_vz.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_vz,
            var_vz.object
        );

    }
    if ((var_dvx.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_dvx,
            var_dvx.object
        );

    }
    if ((var_dvy.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_dvy,
            var_dvy.object
        );

    }
    if ((var_dvz.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_dvz,
            var_dvz.object
        );

    }
    if ((var_twice_max_segment_length.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_twice_max_segment_length,
            var_twice_max_segment_length.object
        );

    }
    if ((var_i.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_i,
            var_i.object
        );

    }
    if ((var_already_jumped.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_already_jumped,
            var_already_jumped.object
        );

    }
    if ((var_point.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_point,
            var_point.object
        );

    }
    if ((var_distance.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_distance,
            var_distance.object
        );

    }
    if ((var_next_point.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_next_point,
            var_next_point.object
        );

    }
    if ((var_jump_x.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_jump_x,
            var_jump_x.object
        );

    }
    if ((var_jump_y.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_jump_y,
            var_jump_y.object
        );

    }
    if ((var_jump_z.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_jump_z,
            var_jump_z.object
        );

    }
    if ((var_intersect.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_intersect,
            var_intersect.object
        );

    }
    if ((var_intersect_i.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_intersect_i,
            var_intersect_i.object
        );

    }
    if ((var_intersect_j.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_intersect_j,
            var_intersect_j.object
        );

    }
    if ((var_pz.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_pz,
            var_pz.object
        );

    }
    if ((var_dpz.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_dpz,
            var_dpz.object
        );

    }
    if ((var_crossing_sign.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_crossing_sign,
            var_crossing_sign.object
        );

    }
    if ((var_crossing_direction.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_crossing_direction,
            var_crossing_direction.object
        );

    }
    if ((var_num_jumps.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_num_jumps,
            var_num_jumps.object
        );

    }
    if ((var_distance_travelled.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_distance_travelled,
            var_distance_travelled.object
        );

    }
    if ((var_jumps.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_jumps,
            var_jumps.object
        );

    }
    if ((par_v.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_v,
            par_v.object
        );

    }
    if ((par_dv.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_dv,
            par_dv.object
        );

    }
    if ((par_points.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_points,
            par_points.object
        );

    }
    if ((par_segment_lengths.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_segment_lengths,
            par_segment_lengths.object
        );

    }
    if ((par_current_index.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_current_index,
            par_current_index.object
        );

    }
    if ((par_comparison_index.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_comparison_index,
            par_comparison_index.object
        );

    }
    if ((par_max_segment_length.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_max_segment_length,
            par_max_segment_length.object
        );

    }
    if ((par_jump_mode.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_jump_mode,
            par_jump_mode.object
        );

    }
    detachFrame( exception_tb, tmp_frame_locals );


    popFrameStack();

#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );

    // Return the error.
    goto function_exception_exit;
    frame_no_exception_1:;

    // Return statement must be present.
    assert(false);
function_exception_exit:
    assert( exception_type );
    PyErr_Restore( exception_type, exception_value, (PyObject *)exception_tb );
    return NULL;
function_return_exit:
    return tmp_return_value;

}
static PyObject *fparse_function_1_find_crossings_of_module_helpers( Nuitka_FunctionObject *self, PyObject **args, Py_ssize_t args_size, PyObject *kw )
{
    assert( kw == NULL || PyDict_Check( kw ) );

    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_size = kw ? PyDict_Size( kw ) : 0;
    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_found = 0;
    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_only_found = 0;
    Py_ssize_t args_given = args_size;
    PyObject *_python_par_v = NULL;
    PyObject *_python_par_dv = NULL;
    PyObject *_python_par_points = NULL;
    PyObject *_python_par_segment_lengths = NULL;
    PyObject *_python_par_current_index = NULL;
    PyObject *_python_par_comparison_index = NULL;
    PyObject *_python_par_max_segment_length = NULL;
    PyObject *_python_par_jump_mode = NULL;
    // Copy given dictionary values to the the respective variables:
    if ( kw_size > 0 )
    {
        Py_ssize_t ppos = 0;
        PyObject *key, *value;

        while( PyDict_Next( kw, &ppos, &key, &value ) )
        {
#if PYTHON_VERSION < 300
            if (unlikely( !PyString_Check( key ) && !PyUnicode_Check( key ) ))
#else
            if (unlikely( !PyUnicode_Check( key ) ))
#endif
            {
                PyErr_Format( PyExc_TypeError, "find_crossings() keywords must be strings" );
                goto error_exit;
            }

            NUITKA_MAY_BE_UNUSED bool found = false;

            Py_INCREF( key );
            Py_INCREF( value );

            // Quick path, could be our value.
            if ( found == false && const_str_plain_v == key )
            {
                assert( _python_par_v == NULL );
                _python_par_v = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_dv == key )
            {
                assert( _python_par_dv == NULL );
                _python_par_dv = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_points == key )
            {
                assert( _python_par_points == NULL );
                _python_par_points = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_segment_lengths == key )
            {
                assert( _python_par_segment_lengths == NULL );
                _python_par_segment_lengths = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_current_index == key )
            {
                assert( _python_par_current_index == NULL );
                _python_par_current_index = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_comparison_index == key )
            {
                assert( _python_par_comparison_index == NULL );
                _python_par_comparison_index = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_max_segment_length == key )
            {
                assert( _python_par_max_segment_length == NULL );
                _python_par_max_segment_length = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_jump_mode == key )
            {
                assert( _python_par_jump_mode == NULL );
                _python_par_jump_mode = value;

                found = true;
                kw_found += 1;
            }

            // Slow path, compare against all parameter names.
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_v, key ) == 1 )
            {
                assert( _python_par_v == NULL );
                _python_par_v = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_dv, key ) == 1 )
            {
                assert( _python_par_dv == NULL );
                _python_par_dv = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_points, key ) == 1 )
            {
                assert( _python_par_points == NULL );
                _python_par_points = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_segment_lengths, key ) == 1 )
            {
                assert( _python_par_segment_lengths == NULL );
                _python_par_segment_lengths = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_current_index, key ) == 1 )
            {
                assert( _python_par_current_index == NULL );
                _python_par_current_index = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_comparison_index, key ) == 1 )
            {
                assert( _python_par_comparison_index == NULL );
                _python_par_comparison_index = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_max_segment_length, key ) == 1 )
            {
                assert( _python_par_max_segment_length == NULL );
                _python_par_max_segment_length = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_jump_mode, key ) == 1 )
            {
                assert( _python_par_jump_mode == NULL );
                _python_par_jump_mode = value;

                found = true;
                kw_found += 1;
            }


            Py_DECREF( key );

            if ( found == false )
            {
               Py_DECREF( value );

               PyErr_Format(
                   PyExc_TypeError,
                   "find_crossings() got an unexpected keyword argument '%s'",
                   Nuitka_String_Check( key ) ? Nuitka_String_AsString( key ) : "<non-string>"
               );

               goto error_exit;
            }
        }

#if PYTHON_VERSION < 300
        assert( kw_found == kw_size );
        assert( kw_only_found == 0 );
#endif
    }

    // Check if too many arguments were given in case of non star args
    if (unlikely( args_given > 8 ))
    {
#if PYTHON_VERSION < 270
        ERROR_TOO_MANY_ARGUMENTS( self, args_given, kw_size );
#elif PYTHON_VERSION < 330
        ERROR_TOO_MANY_ARGUMENTS( self, args_given + kw_found );
#else
        ERROR_TOO_MANY_ARGUMENTS( self, args_given, kw_only_found );
#endif
        goto error_exit;
    }


    // Copy normal parameter values given as part of the args list to the respective variables:

    if (likely( 0 < args_given ))
    {
         if (unlikely( _python_par_v != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 0 );
             goto error_exit;
         }

        _python_par_v = INCREASE_REFCOUNT( args[ 0 ] );
    }
    else if ( _python_par_v == NULL )
    {
        if ( 0 + self->m_defaults_given >= 8  )
        {
            _python_par_v = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 0 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 1 < args_given ))
    {
         if (unlikely( _python_par_dv != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 1 );
             goto error_exit;
         }

        _python_par_dv = INCREASE_REFCOUNT( args[ 1 ] );
    }
    else if ( _python_par_dv == NULL )
    {
        if ( 1 + self->m_defaults_given >= 8  )
        {
            _python_par_dv = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 1 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 2 < args_given ))
    {
         if (unlikely( _python_par_points != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 2 );
             goto error_exit;
         }

        _python_par_points = INCREASE_REFCOUNT( args[ 2 ] );
    }
    else if ( _python_par_points == NULL )
    {
        if ( 2 + self->m_defaults_given >= 8  )
        {
            _python_par_points = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 2 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 3 < args_given ))
    {
         if (unlikely( _python_par_segment_lengths != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 3 );
             goto error_exit;
         }

        _python_par_segment_lengths = INCREASE_REFCOUNT( args[ 3 ] );
    }
    else if ( _python_par_segment_lengths == NULL )
    {
        if ( 3 + self->m_defaults_given >= 8  )
        {
            _python_par_segment_lengths = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 3 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 4 < args_given ))
    {
         if (unlikely( _python_par_current_index != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 4 );
             goto error_exit;
         }

        _python_par_current_index = INCREASE_REFCOUNT( args[ 4 ] );
    }
    else if ( _python_par_current_index == NULL )
    {
        if ( 4 + self->m_defaults_given >= 8  )
        {
            _python_par_current_index = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 4 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 5 < args_given ))
    {
         if (unlikely( _python_par_comparison_index != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 5 );
             goto error_exit;
         }

        _python_par_comparison_index = INCREASE_REFCOUNT( args[ 5 ] );
    }
    else if ( _python_par_comparison_index == NULL )
    {
        if ( 5 + self->m_defaults_given >= 8  )
        {
            _python_par_comparison_index = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 5 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 6 < args_given ))
    {
         if (unlikely( _python_par_max_segment_length != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 6 );
             goto error_exit;
         }

        _python_par_max_segment_length = INCREASE_REFCOUNT( args[ 6 ] );
    }
    else if ( _python_par_max_segment_length == NULL )
    {
        if ( 6 + self->m_defaults_given >= 8  )
        {
            _python_par_max_segment_length = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 6 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 7 < args_given ))
    {
         if (unlikely( _python_par_jump_mode != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 7 );
             goto error_exit;
         }

        _python_par_jump_mode = INCREASE_REFCOUNT( args[ 7 ] );
    }
    else if ( _python_par_jump_mode == NULL )
    {
        if ( 7 + self->m_defaults_given >= 8  )
        {
            _python_par_jump_mode = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 7 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }

#if PYTHON_VERSION >= 330
    if (unlikely( _python_par_v == NULL || _python_par_dv == NULL || _python_par_points == NULL || _python_par_segment_lengths == NULL || _python_par_current_index == NULL || _python_par_comparison_index == NULL || _python_par_max_segment_length == NULL || _python_par_jump_mode == NULL ))
    {
        PyObject *values[] = { _python_par_v, _python_par_dv, _python_par_points, _python_par_segment_lengths, _python_par_current_index, _python_par_comparison_index, _python_par_max_segment_length, _python_par_jump_mode };
        ERROR_TOO_FEW_ARGUMENTS( self, values );

        goto error_exit;
    }
#endif


    return impl_function_1_find_crossings_of_module_helpers( self, _python_par_v, _python_par_dv, _python_par_points, _python_par_segment_lengths, _python_par_current_index, _python_par_comparison_index, _python_par_max_segment_length, _python_par_jump_mode );

error_exit:;

    Py_XDECREF( _python_par_v );
    Py_XDECREF( _python_par_dv );
    Py_XDECREF( _python_par_points );
    Py_XDECREF( _python_par_segment_lengths );
    Py_XDECREF( _python_par_current_index );
    Py_XDECREF( _python_par_comparison_index );
    Py_XDECREF( _python_par_max_segment_length );
    Py_XDECREF( _python_par_jump_mode );

    return NULL;
}

static PyObject *dparse_function_1_find_crossings_of_module_helpers( Nuitka_FunctionObject *self, PyObject **args, int size )
{
    if ( size == 8 )
    {
        return impl_function_1_find_crossings_of_module_helpers( self, INCREASE_REFCOUNT( args[ 0 ] ), INCREASE_REFCOUNT( args[ 1 ] ), INCREASE_REFCOUNT( args[ 2 ] ), INCREASE_REFCOUNT( args[ 3 ] ), INCREASE_REFCOUNT( args[ 4 ] ), INCREASE_REFCOUNT( args[ 5 ] ), INCREASE_REFCOUNT( args[ 6 ] ), INCREASE_REFCOUNT( args[ 7 ] ) );
    }
    else
    {
        PyObject *result = fparse_function_1_find_crossings_of_module_helpers( self, args, size, NULL );
        return result;
    }

}



static PyObject *impl_function_2_do_vectors_intersect_of_module_helpers( Nuitka_FunctionObject *self, PyObject *_python_par_px, PyObject *_python_par_py, PyObject *_python_par_dpx, PyObject *_python_par_dpy, PyObject *_python_par_qx, PyObject *_python_par_qy, PyObject *_python_par_dqx, PyObject *_python_par_dqy )
{
    // No context is used.

    // Local variable declarations.
    PyObjectLocalVariable par_px; par_px.object = _python_par_px;
    PyObjectLocalVariable par_py; par_py.object = _python_par_py;
    PyObjectLocalVariable par_dpx; par_dpx.object = _python_par_dpx;
    PyObjectLocalVariable par_dpy; par_dpy.object = _python_par_dpy;
    PyObjectLocalVariable par_qx; par_qx.object = _python_par_qx;
    PyObjectLocalVariable par_qy; par_qy.object = _python_par_qy;
    PyObjectLocalVariable par_dqx; par_dqx.object = _python_par_dqx;
    PyObjectLocalVariable par_dqy; par_dqy.object = _python_par_dqy;
    PyObjectLocalVariable var_t;
    PyObjectLocalVariable var_u;
    PyObjectTempVariable tmp_and_1__value_1;
    PyObjectTempVariable tmp_and_2__value_1;
    PyObject *exception_type = NULL, *exception_value = NULL;
    PyTracebackObject *exception_tb = NULL;
    PyObject *exception_keeper_type_1;
    PyObject *exception_keeper_value_1;
    PyTracebackObject *exception_keeper_tb_1;
    PyObject *exception_keeper_type_2;
    PyObject *exception_keeper_value_2;
    PyTracebackObject *exception_keeper_tb_2;
    PyObject *exception_keeper_type_3;
    PyObject *exception_keeper_value_3;
    PyTracebackObject *exception_keeper_tb_3;
    PyObject *exception_keeper_type_4;
    PyObject *exception_keeper_value_4;
    PyTracebackObject *exception_keeper_tb_4;
    PyObject *exception_keeper_type_5;
    PyObject *exception_keeper_value_5;
    PyTracebackObject *exception_keeper_tb_5;
    PyObject *exception_keeper_type_6;
    PyObject *exception_keeper_value_6;
    PyTracebackObject *exception_keeper_tb_6;
    PyObject *tmp_assign_source_1;
    PyObject *tmp_assign_source_2;
    PyObject *tmp_assign_source_3;
    PyObject *tmp_assign_source_4;
    PyObject *tmp_binop_left_1;
    PyObject *tmp_binop_left_2;
    PyObject *tmp_binop_left_3;
    PyObject *tmp_binop_left_4;
    PyObject *tmp_binop_left_5;
    PyObject *tmp_binop_left_6;
    PyObject *tmp_binop_right_1;
    PyObject *tmp_binop_right_2;
    PyObject *tmp_binop_right_3;
    PyObject *tmp_binop_right_4;
    PyObject *tmp_binop_right_5;
    PyObject *tmp_binop_right_6;
    PyObject *tmp_call_arg_element_1;
    PyObject *tmp_call_arg_element_2;
    PyObject *tmp_call_arg_element_3;
    PyObject *tmp_call_arg_element_4;
    PyObject *tmp_call_arg_element_5;
    PyObject *tmp_call_arg_element_6;
    PyObject *tmp_call_arg_element_7;
    PyObject *tmp_call_arg_element_8;
    PyObject *tmp_call_arg_element_9;
    PyObject *tmp_call_arg_element_10;
    PyObject *tmp_call_arg_element_11;
    PyObject *tmp_call_arg_element_12;
    PyObject *tmp_call_arg_element_13;
    PyObject *tmp_call_arg_element_14;
    PyObject *tmp_call_arg_element_15;
    PyObject *tmp_call_arg_element_16;
    PyObject *tmp_call_arg_element_17;
    PyObject *tmp_call_arg_element_18;
    PyObject *tmp_call_arg_element_19;
    PyObject *tmp_call_arg_element_20;
    PyObject *tmp_call_arg_element_21;
    PyObject *tmp_called_1;
    PyObject *tmp_called_2;
    PyObject *tmp_called_3;
    PyObject *tmp_called_4;
    PyObject *tmp_called_5;
    PyObject *tmp_called_6;
    int tmp_cmp_Lt_1;
    PyObject *tmp_compare_left_1;
    PyObject *tmp_compare_right_1;
    PyObject *tmp_compexpr_left_1;
    PyObject *tmp_compexpr_left_2;
    PyObject *tmp_compexpr_left_3;
    PyObject *tmp_compexpr_left_4;
    PyObject *tmp_compexpr_right_1;
    PyObject *tmp_compexpr_right_2;
    PyObject *tmp_compexpr_right_3;
    PyObject *tmp_compexpr_right_4;
    int tmp_cond_truth_1;
    int tmp_cond_truth_2;
    int tmp_cond_truth_3;
    int tmp_cond_truth_4;
    PyObject *tmp_cond_value_1;
    PyObject *tmp_cond_value_2;
    PyObject *tmp_cond_value_3;
    PyObject *tmp_cond_value_4;
    PyObject *tmp_frame_locals;
    bool tmp_result;
    PyObject *tmp_return_value;
    int tmp_tried_lineno_1;
    int tmp_tried_lineno_2;
    PyObject *tmp_tuple_element_1;
    tmp_return_value = NULL;

    // Actual function code.
    static PyFrameObject *cache_frame_function = NULL;
    MAKE_OR_REUSE_FRAME( cache_frame_function, codeobj_38efff3ec8c2b4fd2a37e20a63161fd5, module_helpers );
    PyFrameObject *frame_function = cache_frame_function;

    // Push the new frame as the currently active one.
    pushFrameStack( frame_function );

    // Mark the frame object as in use, ref count 1 will be up for reuse.
    Py_INCREF( frame_function );
    assert( Py_REFCNT( frame_function ) == 2 ); // Frame stack

#if PYTHON_VERSION >= 340
    frame_function->f_executing += 1;
#endif

    // Framed code:
    tmp_called_1 = LOOKUP_BUILTIN( const_str_plain_abs );
    if ( tmp_called_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 127;
        goto frame_exception_exit_1;
    }
    tmp_called_2 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_cross_product );

    if (unlikely( tmp_called_2 == NULL ))
    {
        tmp_called_2 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_cross_product );
    }

    if ( tmp_called_2 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1020 ], 42, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 127;
        goto frame_exception_exit_1;
    }

    tmp_call_arg_element_2 = par_dpx.object;

    if ( tmp_call_arg_element_2 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1661 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 127;
        goto frame_exception_exit_1;
    }

    tmp_call_arg_element_3 = par_dpy.object;

    if ( tmp_call_arg_element_3 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1710 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 127;
        goto frame_exception_exit_1;
    }

    tmp_call_arg_element_4 = par_dqx.object;

    if ( tmp_call_arg_element_4 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1759 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 127;
        goto frame_exception_exit_1;
    }

    tmp_call_arg_element_5 = par_dqy.object;

    if ( tmp_call_arg_element_5 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1808 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 127;
        goto frame_exception_exit_1;
    }

    frame_function->f_lineno = 127;
    tmp_call_arg_element_1 = CALL_FUNCTION_WITH_ARGS4( tmp_called_2, tmp_call_arg_element_2, tmp_call_arg_element_3, tmp_call_arg_element_4, tmp_call_arg_element_5 );
    if ( tmp_call_arg_element_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 127;
        goto frame_exception_exit_1;
    }
    frame_function->f_lineno = 127;
    tmp_compare_left_1 = CALL_FUNCTION_WITH_ARGS1( tmp_called_1, tmp_call_arg_element_1 );
    Py_DECREF( tmp_call_arg_element_1 );
    if ( tmp_compare_left_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 127;
        goto frame_exception_exit_1;
    }
    tmp_compare_right_1 = const_float_1e_minus_05;
    tmp_cmp_Lt_1 = RICH_COMPARE_BOOL_LT( tmp_compare_left_1, tmp_compare_right_1 );
    if ( tmp_cmp_Lt_1 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_compare_left_1 );

        frame_function->f_lineno = 127;
        goto frame_exception_exit_1;
    }
    Py_DECREF( tmp_compare_left_1 );
    if (tmp_cmp_Lt_1 == 1)
    {
        goto branch_yes_1;
    }
    else
    {
        goto branch_no_1;
    }
    branch_yes_1:;
    tmp_return_value = const_tuple_int_0_float_0_0_float_0_0_tuple;
    Py_INCREF( tmp_return_value );
    goto frame_return_exit_1;
    branch_no_1:;
    tmp_called_3 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_cross_product );

    if (unlikely( tmp_called_3 == NULL ))
    {
        tmp_called_3 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_cross_product );
    }

    if ( tmp_called_3 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1020 ], 42, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }

    tmp_binop_left_2 = par_qx.object;

    if ( tmp_binop_left_2 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1857 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }

    tmp_binop_right_2 = par_px.object;

    if ( tmp_binop_right_2 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1905 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }

    tmp_call_arg_element_6 = BINARY_OPERATION_SUB( tmp_binop_left_2, tmp_binop_right_2 );
    if ( tmp_call_arg_element_6 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }
    tmp_binop_left_3 = par_qy.object;

    if ( tmp_binop_left_3 == NULL )
    {
        Py_DECREF( tmp_call_arg_element_6 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1953 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }

    tmp_binop_right_3 = par_py.object;

    if ( tmp_binop_right_3 == NULL )
    {
        Py_DECREF( tmp_call_arg_element_6 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 2001 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }

    tmp_call_arg_element_7 = BINARY_OPERATION_SUB( tmp_binop_left_3, tmp_binop_right_3 );
    if ( tmp_call_arg_element_7 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_call_arg_element_6 );

        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }
    tmp_call_arg_element_8 = par_dqx.object;

    if ( tmp_call_arg_element_8 == NULL )
    {
        Py_DECREF( tmp_call_arg_element_6 );
        Py_DECREF( tmp_call_arg_element_7 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1759 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }

    tmp_call_arg_element_9 = par_dqy.object;

    if ( tmp_call_arg_element_9 == NULL )
    {
        Py_DECREF( tmp_call_arg_element_6 );
        Py_DECREF( tmp_call_arg_element_7 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1808 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }

    frame_function->f_lineno = 130;
    tmp_binop_left_1 = CALL_FUNCTION_WITH_ARGS4( tmp_called_3, tmp_call_arg_element_6, tmp_call_arg_element_7, tmp_call_arg_element_8, tmp_call_arg_element_9 );
    Py_DECREF( tmp_call_arg_element_6 );
    Py_DECREF( tmp_call_arg_element_7 );
    if ( tmp_binop_left_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }
    tmp_called_4 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_cross_product );

    if (unlikely( tmp_called_4 == NULL ))
    {
        tmp_called_4 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_cross_product );
    }

    if ( tmp_called_4 == NULL )
    {
        Py_DECREF( tmp_binop_left_1 );
        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1020 ], 42, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }

    tmp_call_arg_element_10 = par_dpx.object;

    if ( tmp_call_arg_element_10 == NULL )
    {
        Py_DECREF( tmp_binop_left_1 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1661 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }

    tmp_call_arg_element_11 = par_dpy.object;

    if ( tmp_call_arg_element_11 == NULL )
    {
        Py_DECREF( tmp_binop_left_1 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1710 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }

    tmp_call_arg_element_12 = par_dqx.object;

    if ( tmp_call_arg_element_12 == NULL )
    {
        Py_DECREF( tmp_binop_left_1 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1759 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 131;
        goto frame_exception_exit_1;
    }

    tmp_call_arg_element_13 = par_dqy.object;

    if ( tmp_call_arg_element_13 == NULL )
    {
        Py_DECREF( tmp_binop_left_1 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1808 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 131;
        goto frame_exception_exit_1;
    }

    frame_function->f_lineno = 131;
    tmp_binop_right_1 = CALL_FUNCTION_WITH_ARGS4( tmp_called_4, tmp_call_arg_element_10, tmp_call_arg_element_11, tmp_call_arg_element_12, tmp_call_arg_element_13 );
    if ( tmp_binop_right_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_1 );

        frame_function->f_lineno = 131;
        goto frame_exception_exit_1;
    }
    tmp_assign_source_1 = BINARY_OPERATION_DIV( tmp_binop_left_1, tmp_binop_right_1 );
    Py_DECREF( tmp_binop_left_1 );
    Py_DECREF( tmp_binop_right_1 );
    if ( tmp_assign_source_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 130;
        goto frame_exception_exit_1;
    }
    assert( var_t.object == NULL );
    var_t.object = tmp_assign_source_1;

    // Tried code
    tmp_cond_value_1 = NULL;
    // Tried code
    tmp_compexpr_left_1 = var_t.object;

    tmp_compexpr_right_1 = const_float_1_0;
    tmp_assign_source_2 = RICH_COMPARE_LT( tmp_compexpr_left_1, tmp_compexpr_right_1 );
    if ( tmp_assign_source_2 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 132;
        goto try_finally_handler_2;
    }
    assert( tmp_and_1__value_1.object == NULL );
    tmp_and_1__value_1.object = tmp_assign_source_2;

    tmp_cond_value_2 = tmp_and_1__value_1.object;

    tmp_cond_truth_2 = CHECK_IF_TRUE( tmp_cond_value_2 );
    if ( tmp_cond_truth_2 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 132;
        goto try_finally_handler_2;
    }
    if (tmp_cond_truth_2 == 1)
    {
        goto condexpr_true_1;
    }
    else
    {
        goto condexpr_false_1;
    }
    condexpr_true_1:;
    tmp_cond_value_1 = NULL;
    // Tried code
    tmp_result = tmp_and_1__value_1.object != NULL;
    if ( tmp_result == true )
    {
        Py_DECREF( tmp_and_1__value_1.object );
        tmp_and_1__value_1.object = NULL;
    }

    assert( tmp_result != false );
    tmp_compexpr_left_2 = var_t.object;

    tmp_compexpr_right_2 = const_float_0_0;
    tmp_cond_value_1 = RICH_COMPARE_GT( tmp_compexpr_left_2, tmp_compexpr_right_2 );
    if ( tmp_cond_value_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 132;
        goto try_finally_handler_3;
    }
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_3:;
    exception_keeper_type_1 = exception_type;
    exception_keeper_value_1 = exception_value;
    exception_keeper_tb_1 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_1 != NULL )
    {
        exception_type = exception_keeper_type_1;
        exception_value = exception_keeper_value_1;
        exception_tb = exception_keeper_tb_1;

        goto try_finally_handler_2;
    }

    goto finally_end_1;
    finally_end_1:;
    goto condexpr_end_1;
    condexpr_false_1:;
    tmp_cond_value_1 = tmp_and_1__value_1.object;

    Py_INCREF( tmp_cond_value_1 );
    condexpr_end_1:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_2:;
    exception_keeper_type_2 = exception_type;
    exception_keeper_value_2 = exception_value;
    exception_keeper_tb_2 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_2 != NULL )
    {
        exception_type = exception_keeper_type_2;
        exception_value = exception_keeper_value_2;
        exception_tb = exception_keeper_tb_2;

        goto try_finally_handler_1;
    }

    goto finally_end_2;
    finally_end_2:;
    tmp_cond_truth_1 = CHECK_IF_TRUE( tmp_cond_value_1 );
    if ( tmp_cond_truth_1 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_cond_value_1 );

        frame_function->f_lineno = 132;
        goto try_finally_handler_1;
    }
    Py_DECREF( tmp_cond_value_1 );
    if (tmp_cond_truth_1 == 1)
    {
        goto branch_yes_2;
    }
    else
    {
        goto branch_no_2;
    }
    branch_yes_2:;
    tmp_called_5 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_cross_product );

    if (unlikely( tmp_called_5 == NULL ))
    {
        tmp_called_5 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_cross_product );
    }

    if ( tmp_called_5 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1020 ], 42, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }

    tmp_binop_left_5 = par_qx.object;

    if ( tmp_binop_left_5 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1857 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }

    tmp_binop_right_5 = par_px.object;

    if ( tmp_binop_right_5 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1905 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }

    tmp_call_arg_element_14 = BINARY_OPERATION_SUB( tmp_binop_left_5, tmp_binop_right_5 );
    if ( tmp_call_arg_element_14 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }
    tmp_binop_left_6 = par_qy.object;

    if ( tmp_binop_left_6 == NULL )
    {
        Py_DECREF( tmp_call_arg_element_14 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1953 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }

    tmp_binop_right_6 = par_py.object;

    if ( tmp_binop_right_6 == NULL )
    {
        Py_DECREF( tmp_call_arg_element_14 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 2001 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }

    tmp_call_arg_element_15 = BINARY_OPERATION_SUB( tmp_binop_left_6, tmp_binop_right_6 );
    if ( tmp_call_arg_element_15 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_call_arg_element_14 );

        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }
    tmp_call_arg_element_16 = par_dpx.object;

    if ( tmp_call_arg_element_16 == NULL )
    {
        Py_DECREF( tmp_call_arg_element_14 );
        Py_DECREF( tmp_call_arg_element_15 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1661 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }

    tmp_call_arg_element_17 = par_dpy.object;

    if ( tmp_call_arg_element_17 == NULL )
    {
        Py_DECREF( tmp_call_arg_element_14 );
        Py_DECREF( tmp_call_arg_element_15 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1710 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }

    frame_function->f_lineno = 133;
    tmp_binop_left_4 = CALL_FUNCTION_WITH_ARGS4( tmp_called_5, tmp_call_arg_element_14, tmp_call_arg_element_15, tmp_call_arg_element_16, tmp_call_arg_element_17 );
    Py_DECREF( tmp_call_arg_element_14 );
    Py_DECREF( tmp_call_arg_element_15 );
    if ( tmp_binop_left_4 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }
    tmp_called_6 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_cross_product );

    if (unlikely( tmp_called_6 == NULL ))
    {
        tmp_called_6 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_cross_product );
    }

    if ( tmp_called_6 == NULL )
    {
        Py_DECREF( tmp_binop_left_4 );
        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1020 ], 42, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }

    tmp_call_arg_element_18 = par_dpx.object;

    if ( tmp_call_arg_element_18 == NULL )
    {
        Py_DECREF( tmp_binop_left_4 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1661 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }

    tmp_call_arg_element_19 = par_dpy.object;

    if ( tmp_call_arg_element_19 == NULL )
    {
        Py_DECREF( tmp_binop_left_4 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1710 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }

    tmp_call_arg_element_20 = par_dqx.object;

    if ( tmp_call_arg_element_20 == NULL )
    {
        Py_DECREF( tmp_binop_left_4 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1759 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 134;
        goto try_finally_handler_1;
    }

    tmp_call_arg_element_21 = par_dqy.object;

    if ( tmp_call_arg_element_21 == NULL )
    {
        Py_DECREF( tmp_binop_left_4 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1808 ], 49, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 134;
        goto try_finally_handler_1;
    }

    frame_function->f_lineno = 134;
    tmp_binop_right_4 = CALL_FUNCTION_WITH_ARGS4( tmp_called_6, tmp_call_arg_element_18, tmp_call_arg_element_19, tmp_call_arg_element_20, tmp_call_arg_element_21 );
    if ( tmp_binop_right_4 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_4 );

        frame_function->f_lineno = 134;
        goto try_finally_handler_1;
    }
    tmp_assign_source_3 = BINARY_OPERATION_DIV( tmp_binop_left_4, tmp_binop_right_4 );
    Py_DECREF( tmp_binop_left_4 );
    Py_DECREF( tmp_binop_right_4 );
    if ( tmp_assign_source_3 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 133;
        goto try_finally_handler_1;
    }
    assert( var_u.object == NULL );
    var_u.object = tmp_assign_source_3;

    // Tried code
    tmp_cond_value_3 = NULL;
    // Tried code
    tmp_compexpr_left_3 = var_u.object;

    tmp_compexpr_right_3 = const_float_1_0;
    tmp_assign_source_4 = RICH_COMPARE_LT( tmp_compexpr_left_3, tmp_compexpr_right_3 );
    if ( tmp_assign_source_4 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 135;
        goto try_finally_handler_5;
    }
    assert( tmp_and_2__value_1.object == NULL );
    tmp_and_2__value_1.object = tmp_assign_source_4;

    tmp_cond_value_4 = tmp_and_2__value_1.object;

    tmp_cond_truth_4 = CHECK_IF_TRUE( tmp_cond_value_4 );
    if ( tmp_cond_truth_4 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 135;
        goto try_finally_handler_5;
    }
    if (tmp_cond_truth_4 == 1)
    {
        goto condexpr_true_2;
    }
    else
    {
        goto condexpr_false_2;
    }
    condexpr_true_2:;
    tmp_cond_value_3 = NULL;
    // Tried code
    tmp_result = tmp_and_2__value_1.object != NULL;
    if ( tmp_result == true )
    {
        Py_DECREF( tmp_and_2__value_1.object );
        tmp_and_2__value_1.object = NULL;
    }

    assert( tmp_result != false );
    tmp_compexpr_left_4 = var_u.object;

    tmp_compexpr_right_4 = const_float_0_0;
    tmp_cond_value_3 = RICH_COMPARE_GT( tmp_compexpr_left_4, tmp_compexpr_right_4 );
    if ( tmp_cond_value_3 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 135;
        goto try_finally_handler_6;
    }
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_6:;
    exception_keeper_type_3 = exception_type;
    exception_keeper_value_3 = exception_value;
    exception_keeper_tb_3 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_3 != NULL )
    {
        exception_type = exception_keeper_type_3;
        exception_value = exception_keeper_value_3;
        exception_tb = exception_keeper_tb_3;

        goto try_finally_handler_5;
    }

    goto finally_end_3;
    finally_end_3:;
    goto condexpr_end_2;
    condexpr_false_2:;
    tmp_cond_value_3 = tmp_and_2__value_1.object;

    Py_INCREF( tmp_cond_value_3 );
    condexpr_end_2:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_5:;
    exception_keeper_type_4 = exception_type;
    exception_keeper_value_4 = exception_value;
    exception_keeper_tb_4 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_4 != NULL )
    {
        exception_type = exception_keeper_type_4;
        exception_value = exception_keeper_value_4;
        exception_tb = exception_keeper_tb_4;

        goto try_finally_handler_4;
    }

    goto finally_end_4;
    finally_end_4:;
    tmp_cond_truth_3 = CHECK_IF_TRUE( tmp_cond_value_3 );
    if ( tmp_cond_truth_3 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_cond_value_3 );

        frame_function->f_lineno = 135;
        goto try_finally_handler_4;
    }
    Py_DECREF( tmp_cond_value_3 );
    if (tmp_cond_truth_3 == 1)
    {
        goto branch_yes_3;
    }
    else
    {
        goto branch_no_3;
    }
    branch_yes_3:;
    tmp_return_value = PyTuple_New( 3 );
    tmp_tuple_element_1 = const_int_pos_1;
    Py_INCREF( tmp_tuple_element_1 );
    PyTuple_SET_ITEM( tmp_return_value, 0, tmp_tuple_element_1 );
    tmp_tuple_element_1 = var_t.object;

    if ( tmp_tuple_element_1 == NULL )
    {
        Py_DECREF( tmp_return_value );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 2049 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 136;
        goto try_finally_handler_4;
    }

    Py_INCREF( tmp_tuple_element_1 );
    PyTuple_SET_ITEM( tmp_return_value, 1, tmp_tuple_element_1 );
    tmp_tuple_element_1 = var_u.object;

    if ( tmp_tuple_element_1 == NULL )
    {
        Py_DECREF( tmp_return_value );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 2096 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 136;
        goto try_finally_handler_4;
    }

    Py_INCREF( tmp_tuple_element_1 );
    PyTuple_SET_ITEM( tmp_return_value, 2, tmp_tuple_element_1 );
    goto try_finally_handler_start_2;
    branch_no_3:;
    try_finally_handler_start_2:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_4:;
    exception_keeper_type_5 = exception_type;
    exception_keeper_value_5 = exception_value;
    exception_keeper_tb_5 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    tmp_tried_lineno_1 = frame_function->f_lineno;
    Py_XDECREF( tmp_and_2__value_1.object );
    tmp_and_2__value_1.object = NULL;

    frame_function->f_lineno = tmp_tried_lineno_1;
    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_5 != NULL )
    {
        exception_type = exception_keeper_type_5;
        exception_value = exception_keeper_value_5;
        exception_tb = exception_keeper_tb_5;

        goto try_finally_handler_1;
    }

    // Return value if any.
    if ( tmp_return_value != NULL )
    {
        goto try_finally_handler_start_1;
    }

    goto finally_end_5;
    finally_end_5:;
    branch_no_2:;
    try_finally_handler_start_1:;
    // Final block of try/finally
    // Tried block ends with no exception occured, note that.
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    try_finally_handler_1:;
    exception_keeper_type_6 = exception_type;
    exception_keeper_value_6 = exception_value;
    exception_keeper_tb_6 = exception_tb;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;

    tmp_tried_lineno_2 = frame_function->f_lineno;
    Py_XDECREF( tmp_and_1__value_1.object );
    tmp_and_1__value_1.object = NULL;

    frame_function->f_lineno = tmp_tried_lineno_2;
    // Re-reraise as necessary after finally was executed.
    // Reraise exception if any.
    if ( exception_keeper_type_6 != NULL )
    {
        exception_type = exception_keeper_type_6;
        exception_value = exception_keeper_value_6;
        exception_tb = exception_keeper_tb_6;

        goto frame_exception_exit_1;
    }

    // Return value if any.
    if ( tmp_return_value != NULL )
    {
        goto frame_return_exit_1;
    }

    goto finally_end_6;
    finally_end_6:;

#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif
    // Put the previous frame back on top.
    popFrameStack();
#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );
    goto frame_no_exception_1;
    frame_return_exit_1:;
#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif
    popFrameStack();
#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );
    goto function_return_exit;
    frame_exception_exit_1:;
#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif

    if ( exception_tb == NULL )
    {
        exception_tb = MAKE_TRACEBACK( INCREASE_REFCOUNT( frame_function ) );
    }
    else if ( exception_tb->tb_frame != frame_function )
    {
        PyTracebackObject *traceback_new = (PyTracebackObject *)MAKE_TRACEBACK( INCREASE_REFCOUNT( frame_function ) );
        traceback_new->tb_next = exception_tb;
        exception_tb = traceback_new;
    }


    tmp_frame_locals = PyDict_New();
    if ((var_t.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_t,
            var_t.object
        );

    }
    if ((var_u.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_u,
            var_u.object
        );

    }
    if ((par_px.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_px,
            par_px.object
        );

    }
    if ((par_py.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_py,
            par_py.object
        );

    }
    if ((par_dpx.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_dpx,
            par_dpx.object
        );

    }
    if ((par_dpy.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_dpy,
            par_dpy.object
        );

    }
    if ((par_qx.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_qx,
            par_qx.object
        );

    }
    if ((par_qy.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_qy,
            par_qy.object
        );

    }
    if ((par_dqx.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_dqx,
            par_dqx.object
        );

    }
    if ((par_dqy.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_dqy,
            par_dqy.object
        );

    }
    detachFrame( exception_tb, tmp_frame_locals );


    popFrameStack();

#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );

    // Return the error.
    goto function_exception_exit;
    frame_no_exception_1:;
    tmp_return_value = const_tuple_int_0_float__minus_1_0_float__minus_1_0_tuple;
    Py_INCREF( tmp_return_value );
    goto function_return_exit;

    // Return statement must be present.
    assert(false);
function_exception_exit:
    assert( exception_type );
    PyErr_Restore( exception_type, exception_value, (PyObject *)exception_tb );
    return NULL;
function_return_exit:
    return tmp_return_value;

}
static PyObject *fparse_function_2_do_vectors_intersect_of_module_helpers( Nuitka_FunctionObject *self, PyObject **args, Py_ssize_t args_size, PyObject *kw )
{
    assert( kw == NULL || PyDict_Check( kw ) );

    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_size = kw ? PyDict_Size( kw ) : 0;
    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_found = 0;
    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_only_found = 0;
    Py_ssize_t args_given = args_size;
    PyObject *_python_par_px = NULL;
    PyObject *_python_par_py = NULL;
    PyObject *_python_par_dpx = NULL;
    PyObject *_python_par_dpy = NULL;
    PyObject *_python_par_qx = NULL;
    PyObject *_python_par_qy = NULL;
    PyObject *_python_par_dqx = NULL;
    PyObject *_python_par_dqy = NULL;
    // Copy given dictionary values to the the respective variables:
    if ( kw_size > 0 )
    {
        Py_ssize_t ppos = 0;
        PyObject *key, *value;

        while( PyDict_Next( kw, &ppos, &key, &value ) )
        {
#if PYTHON_VERSION < 300
            if (unlikely( !PyString_Check( key ) && !PyUnicode_Check( key ) ))
#else
            if (unlikely( !PyUnicode_Check( key ) ))
#endif
            {
                PyErr_Format( PyExc_TypeError, "do_vectors_intersect() keywords must be strings" );
                goto error_exit;
            }

            NUITKA_MAY_BE_UNUSED bool found = false;

            Py_INCREF( key );
            Py_INCREF( value );

            // Quick path, could be our value.
            if ( found == false && const_str_plain_px == key )
            {
                assert( _python_par_px == NULL );
                _python_par_px = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_py == key )
            {
                assert( _python_par_py == NULL );
                _python_par_py = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_dpx == key )
            {
                assert( _python_par_dpx == NULL );
                _python_par_dpx = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_dpy == key )
            {
                assert( _python_par_dpy == NULL );
                _python_par_dpy = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_qx == key )
            {
                assert( _python_par_qx == NULL );
                _python_par_qx = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_qy == key )
            {
                assert( _python_par_qy == NULL );
                _python_par_qy = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_dqx == key )
            {
                assert( _python_par_dqx == NULL );
                _python_par_dqx = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_dqy == key )
            {
                assert( _python_par_dqy == NULL );
                _python_par_dqy = value;

                found = true;
                kw_found += 1;
            }

            // Slow path, compare against all parameter names.
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_px, key ) == 1 )
            {
                assert( _python_par_px == NULL );
                _python_par_px = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_py, key ) == 1 )
            {
                assert( _python_par_py == NULL );
                _python_par_py = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_dpx, key ) == 1 )
            {
                assert( _python_par_dpx == NULL );
                _python_par_dpx = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_dpy, key ) == 1 )
            {
                assert( _python_par_dpy == NULL );
                _python_par_dpy = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_qx, key ) == 1 )
            {
                assert( _python_par_qx == NULL );
                _python_par_qx = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_qy, key ) == 1 )
            {
                assert( _python_par_qy == NULL );
                _python_par_qy = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_dqx, key ) == 1 )
            {
                assert( _python_par_dqx == NULL );
                _python_par_dqx = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_dqy, key ) == 1 )
            {
                assert( _python_par_dqy == NULL );
                _python_par_dqy = value;

                found = true;
                kw_found += 1;
            }


            Py_DECREF( key );

            if ( found == false )
            {
               Py_DECREF( value );

               PyErr_Format(
                   PyExc_TypeError,
                   "do_vectors_intersect() got an unexpected keyword argument '%s'",
                   Nuitka_String_Check( key ) ? Nuitka_String_AsString( key ) : "<non-string>"
               );

               goto error_exit;
            }
        }

#if PYTHON_VERSION < 300
        assert( kw_found == kw_size );
        assert( kw_only_found == 0 );
#endif
    }

    // Check if too many arguments were given in case of non star args
    if (unlikely( args_given > 8 ))
    {
#if PYTHON_VERSION < 270
        ERROR_TOO_MANY_ARGUMENTS( self, args_given, kw_size );
#elif PYTHON_VERSION < 330
        ERROR_TOO_MANY_ARGUMENTS( self, args_given + kw_found );
#else
        ERROR_TOO_MANY_ARGUMENTS( self, args_given, kw_only_found );
#endif
        goto error_exit;
    }


    // Copy normal parameter values given as part of the args list to the respective variables:

    if (likely( 0 < args_given ))
    {
         if (unlikely( _python_par_px != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 0 );
             goto error_exit;
         }

        _python_par_px = INCREASE_REFCOUNT( args[ 0 ] );
    }
    else if ( _python_par_px == NULL )
    {
        if ( 0 + self->m_defaults_given >= 8  )
        {
            _python_par_px = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 0 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 1 < args_given ))
    {
         if (unlikely( _python_par_py != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 1 );
             goto error_exit;
         }

        _python_par_py = INCREASE_REFCOUNT( args[ 1 ] );
    }
    else if ( _python_par_py == NULL )
    {
        if ( 1 + self->m_defaults_given >= 8  )
        {
            _python_par_py = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 1 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 2 < args_given ))
    {
         if (unlikely( _python_par_dpx != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 2 );
             goto error_exit;
         }

        _python_par_dpx = INCREASE_REFCOUNT( args[ 2 ] );
    }
    else if ( _python_par_dpx == NULL )
    {
        if ( 2 + self->m_defaults_given >= 8  )
        {
            _python_par_dpx = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 2 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 3 < args_given ))
    {
         if (unlikely( _python_par_dpy != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 3 );
             goto error_exit;
         }

        _python_par_dpy = INCREASE_REFCOUNT( args[ 3 ] );
    }
    else if ( _python_par_dpy == NULL )
    {
        if ( 3 + self->m_defaults_given >= 8  )
        {
            _python_par_dpy = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 3 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 4 < args_given ))
    {
         if (unlikely( _python_par_qx != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 4 );
             goto error_exit;
         }

        _python_par_qx = INCREASE_REFCOUNT( args[ 4 ] );
    }
    else if ( _python_par_qx == NULL )
    {
        if ( 4 + self->m_defaults_given >= 8  )
        {
            _python_par_qx = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 4 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 5 < args_given ))
    {
         if (unlikely( _python_par_qy != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 5 );
             goto error_exit;
         }

        _python_par_qy = INCREASE_REFCOUNT( args[ 5 ] );
    }
    else if ( _python_par_qy == NULL )
    {
        if ( 5 + self->m_defaults_given >= 8  )
        {
            _python_par_qy = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 5 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 6 < args_given ))
    {
         if (unlikely( _python_par_dqx != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 6 );
             goto error_exit;
         }

        _python_par_dqx = INCREASE_REFCOUNT( args[ 6 ] );
    }
    else if ( _python_par_dqx == NULL )
    {
        if ( 6 + self->m_defaults_given >= 8  )
        {
            _python_par_dqx = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 6 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 7 < args_given ))
    {
         if (unlikely( _python_par_dqy != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 7 );
             goto error_exit;
         }

        _python_par_dqy = INCREASE_REFCOUNT( args[ 7 ] );
    }
    else if ( _python_par_dqy == NULL )
    {
        if ( 7 + self->m_defaults_given >= 8  )
        {
            _python_par_dqy = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 7 - 8 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }

#if PYTHON_VERSION >= 330
    if (unlikely( _python_par_px == NULL || _python_par_py == NULL || _python_par_dpx == NULL || _python_par_dpy == NULL || _python_par_qx == NULL || _python_par_qy == NULL || _python_par_dqx == NULL || _python_par_dqy == NULL ))
    {
        PyObject *values[] = { _python_par_px, _python_par_py, _python_par_dpx, _python_par_dpy, _python_par_qx, _python_par_qy, _python_par_dqx, _python_par_dqy };
        ERROR_TOO_FEW_ARGUMENTS( self, values );

        goto error_exit;
    }
#endif


    return impl_function_2_do_vectors_intersect_of_module_helpers( self, _python_par_px, _python_par_py, _python_par_dpx, _python_par_dpy, _python_par_qx, _python_par_qy, _python_par_dqx, _python_par_dqy );

error_exit:;

    Py_XDECREF( _python_par_px );
    Py_XDECREF( _python_par_py );
    Py_XDECREF( _python_par_dpx );
    Py_XDECREF( _python_par_dpy );
    Py_XDECREF( _python_par_qx );
    Py_XDECREF( _python_par_qy );
    Py_XDECREF( _python_par_dqx );
    Py_XDECREF( _python_par_dqy );

    return NULL;
}

static PyObject *dparse_function_2_do_vectors_intersect_of_module_helpers( Nuitka_FunctionObject *self, PyObject **args, int size )
{
    if ( size == 8 )
    {
        return impl_function_2_do_vectors_intersect_of_module_helpers( self, INCREASE_REFCOUNT( args[ 0 ] ), INCREASE_REFCOUNT( args[ 1 ] ), INCREASE_REFCOUNT( args[ 2 ] ), INCREASE_REFCOUNT( args[ 3 ] ), INCREASE_REFCOUNT( args[ 4 ] ), INCREASE_REFCOUNT( args[ 5 ] ), INCREASE_REFCOUNT( args[ 6 ] ), INCREASE_REFCOUNT( args[ 7 ] ) );
    }
    else
    {
        PyObject *result = fparse_function_2_do_vectors_intersect_of_module_helpers( self, args, size, NULL );
        return result;
    }

}



static PyObject *impl_function_3_cross_product_of_module_helpers( Nuitka_FunctionObject *self, PyObject *_python_par_px, PyObject *_python_par_py, PyObject *_python_par_qx, PyObject *_python_par_qy )
{
    // No context is used.

    // Local variable declarations.
    PyObjectLocalVariable par_px; par_px.object = _python_par_px;
    PyObjectLocalVariable par_py; par_py.object = _python_par_py;
    PyObjectLocalVariable par_qx; par_qx.object = _python_par_qx;
    PyObjectLocalVariable par_qy; par_qy.object = _python_par_qy;
    PyObject *exception_type = NULL, *exception_value = NULL;
    PyTracebackObject *exception_tb = NULL;
    PyObject *tmp_binop_left_1;
    PyObject *tmp_binop_left_2;
    PyObject *tmp_binop_left_3;
    PyObject *tmp_binop_right_1;
    PyObject *tmp_binop_right_2;
    PyObject *tmp_binop_right_3;
    PyObject *tmp_frame_locals;
    PyObject *tmp_return_value;
    tmp_return_value = NULL;

    // Actual function code.
    static PyFrameObject *cache_frame_function = NULL;
    MAKE_OR_REUSE_FRAME( cache_frame_function, codeobj_641b42b44ac395c476fb706d09df8c7a, module_helpers );
    PyFrameObject *frame_function = cache_frame_function;

    // Push the new frame as the currently active one.
    pushFrameStack( frame_function );

    // Mark the frame object as in use, ref count 1 will be up for reuse.
    Py_INCREF( frame_function );
    assert( Py_REFCNT( frame_function ) == 2 ); // Frame stack

#if PYTHON_VERSION >= 340
    frame_function->f_executing += 1;
#endif

    // Framed code:
    tmp_binop_left_2 = par_px.object;

    if ( tmp_binop_left_2 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1905 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 142;
        goto frame_exception_exit_1;
    }

    tmp_binop_right_2 = par_qy.object;

    if ( tmp_binop_right_2 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1953 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 142;
        goto frame_exception_exit_1;
    }

    tmp_binop_left_1 = BINARY_OPERATION_MUL( tmp_binop_left_2, tmp_binop_right_2 );
    if ( tmp_binop_left_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 142;
        goto frame_exception_exit_1;
    }
    tmp_binop_left_3 = par_py.object;

    if ( tmp_binop_left_3 == NULL )
    {
        Py_DECREF( tmp_binop_left_1 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 2001 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 142;
        goto frame_exception_exit_1;
    }

    tmp_binop_right_3 = par_qx.object;

    if ( tmp_binop_right_3 == NULL )
    {
        Py_DECREF( tmp_binop_left_1 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 1857 ], 48, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 142;
        goto frame_exception_exit_1;
    }

    tmp_binop_right_1 = BINARY_OPERATION_MUL( tmp_binop_left_3, tmp_binop_right_3 );
    if ( tmp_binop_right_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_1 );

        frame_function->f_lineno = 142;
        goto frame_exception_exit_1;
    }
    tmp_return_value = BINARY_OPERATION_SUB( tmp_binop_left_1, tmp_binop_right_1 );
    Py_DECREF( tmp_binop_left_1 );
    Py_DECREF( tmp_binop_right_1 );
    if ( tmp_return_value == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 142;
        goto frame_exception_exit_1;
    }
    goto frame_return_exit_1;

#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif
    // Put the previous frame back on top.
    popFrameStack();
#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );
    goto frame_no_exception_1;
    frame_return_exit_1:;
#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif
    popFrameStack();
#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );
    goto function_return_exit;
    frame_exception_exit_1:;
#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif

    if ( exception_tb == NULL )
    {
        exception_tb = MAKE_TRACEBACK( INCREASE_REFCOUNT( frame_function ) );
    }
    else if ( exception_tb->tb_frame != frame_function )
    {
        PyTracebackObject *traceback_new = (PyTracebackObject *)MAKE_TRACEBACK( INCREASE_REFCOUNT( frame_function ) );
        traceback_new->tb_next = exception_tb;
        exception_tb = traceback_new;
    }


    tmp_frame_locals = PyDict_New();
    if ((par_px.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_px,
            par_px.object
        );

    }
    if ((par_py.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_py,
            par_py.object
        );

    }
    if ((par_qx.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_qx,
            par_qx.object
        );

    }
    if ((par_qy.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_qy,
            par_qy.object
        );

    }
    detachFrame( exception_tb, tmp_frame_locals );


    popFrameStack();

#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );

    // Return the error.
    goto function_exception_exit;
    frame_no_exception_1:;


    // Return statement must be present.
    assert(false);
function_exception_exit:
    assert( exception_type );
    PyErr_Restore( exception_type, exception_value, (PyObject *)exception_tb );
    return NULL;
function_return_exit:
    return tmp_return_value;

}
static PyObject *fparse_function_3_cross_product_of_module_helpers( Nuitka_FunctionObject *self, PyObject **args, Py_ssize_t args_size, PyObject *kw )
{
    assert( kw == NULL || PyDict_Check( kw ) );

    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_size = kw ? PyDict_Size( kw ) : 0;
    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_found = 0;
    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_only_found = 0;
    Py_ssize_t args_given = args_size;
    PyObject *_python_par_px = NULL;
    PyObject *_python_par_py = NULL;
    PyObject *_python_par_qx = NULL;
    PyObject *_python_par_qy = NULL;
    // Copy given dictionary values to the the respective variables:
    if ( kw_size > 0 )
    {
        Py_ssize_t ppos = 0;
        PyObject *key, *value;

        while( PyDict_Next( kw, &ppos, &key, &value ) )
        {
#if PYTHON_VERSION < 300
            if (unlikely( !PyString_Check( key ) && !PyUnicode_Check( key ) ))
#else
            if (unlikely( !PyUnicode_Check( key ) ))
#endif
            {
                PyErr_Format( PyExc_TypeError, "cross_product() keywords must be strings" );
                goto error_exit;
            }

            NUITKA_MAY_BE_UNUSED bool found = false;

            Py_INCREF( key );
            Py_INCREF( value );

            // Quick path, could be our value.
            if ( found == false && const_str_plain_px == key )
            {
                assert( _python_par_px == NULL );
                _python_par_px = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_py == key )
            {
                assert( _python_par_py == NULL );
                _python_par_py = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_qx == key )
            {
                assert( _python_par_qx == NULL );
                _python_par_qx = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_qy == key )
            {
                assert( _python_par_qy == NULL );
                _python_par_qy = value;

                found = true;
                kw_found += 1;
            }

            // Slow path, compare against all parameter names.
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_px, key ) == 1 )
            {
                assert( _python_par_px == NULL );
                _python_par_px = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_py, key ) == 1 )
            {
                assert( _python_par_py == NULL );
                _python_par_py = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_qx, key ) == 1 )
            {
                assert( _python_par_qx == NULL );
                _python_par_qx = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_qy, key ) == 1 )
            {
                assert( _python_par_qy == NULL );
                _python_par_qy = value;

                found = true;
                kw_found += 1;
            }


            Py_DECREF( key );

            if ( found == false )
            {
               Py_DECREF( value );

               PyErr_Format(
                   PyExc_TypeError,
                   "cross_product() got an unexpected keyword argument '%s'",
                   Nuitka_String_Check( key ) ? Nuitka_String_AsString( key ) : "<non-string>"
               );

               goto error_exit;
            }
        }

#if PYTHON_VERSION < 300
        assert( kw_found == kw_size );
        assert( kw_only_found == 0 );
#endif
    }

    // Check if too many arguments were given in case of non star args
    if (unlikely( args_given > 4 ))
    {
#if PYTHON_VERSION < 270
        ERROR_TOO_MANY_ARGUMENTS( self, args_given, kw_size );
#elif PYTHON_VERSION < 330
        ERROR_TOO_MANY_ARGUMENTS( self, args_given + kw_found );
#else
        ERROR_TOO_MANY_ARGUMENTS( self, args_given, kw_only_found );
#endif
        goto error_exit;
    }


    // Copy normal parameter values given as part of the args list to the respective variables:

    if (likely( 0 < args_given ))
    {
         if (unlikely( _python_par_px != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 0 );
             goto error_exit;
         }

        _python_par_px = INCREASE_REFCOUNT( args[ 0 ] );
    }
    else if ( _python_par_px == NULL )
    {
        if ( 0 + self->m_defaults_given >= 4  )
        {
            _python_par_px = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 0 - 4 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 1 < args_given ))
    {
         if (unlikely( _python_par_py != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 1 );
             goto error_exit;
         }

        _python_par_py = INCREASE_REFCOUNT( args[ 1 ] );
    }
    else if ( _python_par_py == NULL )
    {
        if ( 1 + self->m_defaults_given >= 4  )
        {
            _python_par_py = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 1 - 4 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 2 < args_given ))
    {
         if (unlikely( _python_par_qx != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 2 );
             goto error_exit;
         }

        _python_par_qx = INCREASE_REFCOUNT( args[ 2 ] );
    }
    else if ( _python_par_qx == NULL )
    {
        if ( 2 + self->m_defaults_given >= 4  )
        {
            _python_par_qx = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 2 - 4 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 3 < args_given ))
    {
         if (unlikely( _python_par_qy != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 3 );
             goto error_exit;
         }

        _python_par_qy = INCREASE_REFCOUNT( args[ 3 ] );
    }
    else if ( _python_par_qy == NULL )
    {
        if ( 3 + self->m_defaults_given >= 4  )
        {
            _python_par_qy = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 3 - 4 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }

#if PYTHON_VERSION >= 330
    if (unlikely( _python_par_px == NULL || _python_par_py == NULL || _python_par_qx == NULL || _python_par_qy == NULL ))
    {
        PyObject *values[] = { _python_par_px, _python_par_py, _python_par_qx, _python_par_qy };
        ERROR_TOO_FEW_ARGUMENTS( self, values );

        goto error_exit;
    }
#endif


    return impl_function_3_cross_product_of_module_helpers( self, _python_par_px, _python_par_py, _python_par_qx, _python_par_qy );

error_exit:;

    Py_XDECREF( _python_par_px );
    Py_XDECREF( _python_par_py );
    Py_XDECREF( _python_par_qx );
    Py_XDECREF( _python_par_qy );

    return NULL;
}

static PyObject *dparse_function_3_cross_product_of_module_helpers( Nuitka_FunctionObject *self, PyObject **args, int size )
{
    if ( size == 4 )
    {
        return impl_function_3_cross_product_of_module_helpers( self, INCREASE_REFCOUNT( args[ 0 ] ), INCREASE_REFCOUNT( args[ 1 ] ), INCREASE_REFCOUNT( args[ 2 ] ), INCREASE_REFCOUNT( args[ 3 ] ) );
    }
    else
    {
        PyObject *result = fparse_function_3_cross_product_of_module_helpers( self, args, size, NULL );
        return result;
    }

}



static PyObject *impl_function_4_sign_of_module_helpers( Nuitka_FunctionObject *self, PyObject *_python_par_a )
{
    // No context is used.

    // Local variable declarations.
    PyObjectLocalVariable par_a; par_a.object = _python_par_a;
    PyObject *exception_type = NULL, *exception_value = NULL;
    PyTracebackObject *exception_tb = NULL;
    int tmp_cmp_Gt_1;
    int tmp_cmp_Lt_1;
    PyObject *tmp_compare_left_1;
    PyObject *tmp_compare_left_2;
    PyObject *tmp_compare_right_1;
    PyObject *tmp_compare_right_2;
    PyObject *tmp_frame_locals;
    PyObject *tmp_return_value;
    tmp_return_value = NULL;

    // Actual function code.
    static PyFrameObject *cache_frame_function = NULL;
    MAKE_OR_REUSE_FRAME( cache_frame_function, codeobj_fd2d6789440b495c075e8467bf2253e1, module_helpers );
    PyFrameObject *frame_function = cache_frame_function;

    // Push the new frame as the currently active one.
    pushFrameStack( frame_function );

    // Mark the frame object as in use, ref count 1 will be up for reuse.
    Py_INCREF( frame_function );
    assert( Py_REFCNT( frame_function ) == 2 ); // Frame stack

#if PYTHON_VERSION >= 340
    frame_function->f_executing += 1;
#endif

    // Framed code:
    tmp_compare_left_1 = par_a.object;

    if ( tmp_compare_left_1 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 2143 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 146;
        goto frame_exception_exit_1;
    }

    tmp_compare_right_1 = const_float_0_0;
    tmp_cmp_Gt_1 = RICH_COMPARE_BOOL_GT( tmp_compare_left_1, tmp_compare_right_1 );
    if ( tmp_cmp_Gt_1 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 146;
        goto frame_exception_exit_1;
    }
    if (tmp_cmp_Gt_1 == 1)
    {
        goto condexpr_true_1;
    }
    else
    {
        goto condexpr_false_1;
    }
    condexpr_true_1:;
    tmp_return_value = const_float_1_0;
    goto condexpr_end_1;
    condexpr_false_1:;
    tmp_compare_left_2 = par_a.object;

    if ( tmp_compare_left_2 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 2143 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 146;
        goto frame_exception_exit_1;
    }

    tmp_compare_right_2 = const_float_0_0;
    tmp_cmp_Lt_1 = RICH_COMPARE_BOOL_LT( tmp_compare_left_2, tmp_compare_right_2 );
    if ( tmp_cmp_Lt_1 == -1 )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 146;
        goto frame_exception_exit_1;
    }
    if (tmp_cmp_Lt_1 == 1)
    {
        goto condexpr_true_2;
    }
    else
    {
        goto condexpr_false_2;
    }
    condexpr_true_2:;
    tmp_return_value = const_float__minus_1_0;
    goto condexpr_end_2;
    condexpr_false_2:;
    tmp_return_value = const_float_0_0;
    condexpr_end_2:;
    condexpr_end_1:;
    Py_INCREF( tmp_return_value );
    goto frame_return_exit_1;

#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif
    // Put the previous frame back on top.
    popFrameStack();
#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );
    goto frame_no_exception_1;
    frame_return_exit_1:;
#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif
    popFrameStack();
#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );
    goto function_return_exit;
    frame_exception_exit_1:;
#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif

    if ( exception_tb == NULL )
    {
        exception_tb = MAKE_TRACEBACK( INCREASE_REFCOUNT( frame_function ) );
    }
    else if ( exception_tb->tb_frame != frame_function )
    {
        PyTracebackObject *traceback_new = (PyTracebackObject *)MAKE_TRACEBACK( INCREASE_REFCOUNT( frame_function ) );
        traceback_new->tb_next = exception_tb;
        exception_tb = traceback_new;
    }


    tmp_frame_locals = PyDict_New();
    if ((par_a.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_a,
            par_a.object
        );

    }
    detachFrame( exception_tb, tmp_frame_locals );


    popFrameStack();

#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );

    // Return the error.
    goto function_exception_exit;
    frame_no_exception_1:;


    // Return statement must be present.
    assert(false);
function_exception_exit:
    assert( exception_type );
    PyErr_Restore( exception_type, exception_value, (PyObject *)exception_tb );
    return NULL;
function_return_exit:
    return tmp_return_value;

}
static PyObject *fparse_function_4_sign_of_module_helpers( Nuitka_FunctionObject *self, PyObject **args, Py_ssize_t args_size, PyObject *kw )
{
    assert( kw == NULL || PyDict_Check( kw ) );

    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_size = kw ? PyDict_Size( kw ) : 0;
    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_found = 0;
    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_only_found = 0;
    Py_ssize_t args_given = args_size;
    PyObject *_python_par_a = NULL;
    // Copy given dictionary values to the the respective variables:
    if ( kw_size > 0 )
    {
        Py_ssize_t ppos = 0;
        PyObject *key, *value;

        while( PyDict_Next( kw, &ppos, &key, &value ) )
        {
#if PYTHON_VERSION < 300
            if (unlikely( !PyString_Check( key ) && !PyUnicode_Check( key ) ))
#else
            if (unlikely( !PyUnicode_Check( key ) ))
#endif
            {
                PyErr_Format( PyExc_TypeError, "sign() keywords must be strings" );
                goto error_exit;
            }

            NUITKA_MAY_BE_UNUSED bool found = false;

            Py_INCREF( key );
            Py_INCREF( value );

            // Quick path, could be our value.
            if ( found == false && const_str_plain_a == key )
            {
                assert( _python_par_a == NULL );
                _python_par_a = value;

                found = true;
                kw_found += 1;
            }

            // Slow path, compare against all parameter names.
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_a, key ) == 1 )
            {
                assert( _python_par_a == NULL );
                _python_par_a = value;

                found = true;
                kw_found += 1;
            }


            Py_DECREF( key );

            if ( found == false )
            {
               Py_DECREF( value );

               PyErr_Format(
                   PyExc_TypeError,
                   "sign() got an unexpected keyword argument '%s'",
                   Nuitka_String_Check( key ) ? Nuitka_String_AsString( key ) : "<non-string>"
               );

               goto error_exit;
            }
        }

#if PYTHON_VERSION < 300
        assert( kw_found == kw_size );
        assert( kw_only_found == 0 );
#endif
    }

    // Check if too many arguments were given in case of non star args
    if (unlikely( args_given > 1 ))
    {
#if PYTHON_VERSION < 270
        ERROR_TOO_MANY_ARGUMENTS( self, args_given, kw_size );
#elif PYTHON_VERSION < 330
        ERROR_TOO_MANY_ARGUMENTS( self, args_given + kw_found );
#else
        ERROR_TOO_MANY_ARGUMENTS( self, args_given, kw_only_found );
#endif
        goto error_exit;
    }


    // Copy normal parameter values given as part of the args list to the respective variables:

    if (likely( 0 < args_given ))
    {
         if (unlikely( _python_par_a != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 0 );
             goto error_exit;
         }

        _python_par_a = INCREASE_REFCOUNT( args[ 0 ] );
    }
    else if ( _python_par_a == NULL )
    {
        if ( 0 + self->m_defaults_given >= 1  )
        {
            _python_par_a = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 0 - 1 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }

#if PYTHON_VERSION >= 330
    if (unlikely( _python_par_a == NULL ))
    {
        PyObject *values[] = { _python_par_a };
        ERROR_TOO_FEW_ARGUMENTS( self, values );

        goto error_exit;
    }
#endif


    return impl_function_4_sign_of_module_helpers( self, _python_par_a );

error_exit:;

    Py_XDECREF( _python_par_a );

    return NULL;
}

static PyObject *dparse_function_4_sign_of_module_helpers( Nuitka_FunctionObject *self, PyObject **args, int size )
{
    if ( size == 1 )
    {
        return impl_function_4_sign_of_module_helpers( self, INCREASE_REFCOUNT( args[ 0 ] ) );
    }
    else
    {
        PyObject *result = fparse_function_4_sign_of_module_helpers( self, args, size, NULL );
        return result;
    }

}



static PyObject *impl_function_5_mag_difference_of_module_helpers( Nuitka_FunctionObject *self, PyObject *_python_par_a, PyObject *_python_par_b )
{
    // No context is used.

    // Local variable declarations.
    PyObjectLocalVariable par_a; par_a.object = _python_par_a;
    PyObjectLocalVariable par_b; par_b.object = _python_par_b;
    PyObject *exception_type = NULL, *exception_value = NULL;
    PyTracebackObject *exception_tb = NULL;
    PyObject *tmp_binop_left_1;
    PyObject *tmp_binop_left_2;
    PyObject *tmp_binop_left_3;
    PyObject *tmp_binop_left_4;
    PyObject *tmp_binop_left_5;
    PyObject *tmp_binop_right_1;
    PyObject *tmp_binop_right_2;
    PyObject *tmp_binop_right_3;
    PyObject *tmp_binop_right_4;
    PyObject *tmp_binop_right_5;
    PyObject *tmp_call_arg_element_1;
    PyObject *tmp_called_1;
    PyObject *tmp_frame_locals;
    PyObject *tmp_return_value;
    PyObject *tmp_subscr_subscript_1;
    PyObject *tmp_subscr_subscript_2;
    PyObject *tmp_subscr_subscript_3;
    PyObject *tmp_subscr_subscript_4;
    PyObject *tmp_subscr_target_1;
    PyObject *tmp_subscr_target_2;
    PyObject *tmp_subscr_target_3;
    PyObject *tmp_subscr_target_4;
    tmp_return_value = NULL;

    // Actual function code.
    static PyFrameObject *cache_frame_function = NULL;
    MAKE_OR_REUSE_FRAME( cache_frame_function, codeobj_71cfbdc35b90f7d090c4273357868645, module_helpers );
    PyFrameObject *frame_function = cache_frame_function;

    // Push the new frame as the currently active one.
    pushFrameStack( frame_function );

    // Mark the frame object as in use, ref count 1 will be up for reuse.
    Py_INCREF( frame_function );
    assert( Py_REFCNT( frame_function ) == 2 ); // Frame stack

#if PYTHON_VERSION >= 340
    frame_function->f_executing += 1;
#endif

    // Framed code:
    tmp_called_1 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_csqrt );

    if (unlikely( tmp_called_1 == NULL ))
    {
        tmp_called_1 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_csqrt );
    }

    if ( tmp_called_1 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 281 ], 34, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }

    tmp_subscr_target_1 = par_b.object;

    if ( tmp_subscr_target_1 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 2190 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }

    tmp_subscr_subscript_1 = const_int_0;
    tmp_binop_left_3 = LOOKUP_SUBSCRIPT( tmp_subscr_target_1, tmp_subscr_subscript_1 );
    if ( tmp_binop_left_3 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }
    tmp_subscr_target_2 = par_a.object;

    if ( tmp_subscr_target_2 == NULL )
    {
        Py_DECREF( tmp_binop_left_3 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 2143 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }

    tmp_subscr_subscript_2 = const_int_0;
    tmp_binop_right_3 = LOOKUP_SUBSCRIPT( tmp_subscr_target_2, tmp_subscr_subscript_2 );
    if ( tmp_binop_right_3 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_3 );

        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }
    tmp_binop_left_2 = BINARY_OPERATION_SUB( tmp_binop_left_3, tmp_binop_right_3 );
    Py_DECREF( tmp_binop_left_3 );
    Py_DECREF( tmp_binop_right_3 );
    if ( tmp_binop_left_2 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }
    tmp_binop_right_2 = const_int_pos_2;
    tmp_binop_left_1 = POWER_OPERATION( tmp_binop_left_2, tmp_binop_right_2 );
    Py_DECREF( tmp_binop_left_2 );
    if ( tmp_binop_left_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }
    tmp_subscr_target_3 = par_b.object;

    if ( tmp_subscr_target_3 == NULL )
    {
        Py_DECREF( tmp_binop_left_1 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 2190 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }

    tmp_subscr_subscript_3 = const_int_pos_1;
    tmp_binop_left_5 = LOOKUP_SUBSCRIPT( tmp_subscr_target_3, tmp_subscr_subscript_3 );
    if ( tmp_binop_left_5 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_1 );

        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }
    tmp_subscr_target_4 = par_a.object;

    if ( tmp_subscr_target_4 == NULL )
    {
        Py_DECREF( tmp_binop_left_1 );
        Py_DECREF( tmp_binop_left_5 );
        exception_type = INCREASE_REFCOUNT( PyExc_UnboundLocalError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 2143 ], 47, 0 );
        exception_tb = NULL;

        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }

    tmp_subscr_subscript_4 = const_int_pos_1;
    tmp_binop_right_5 = LOOKUP_SUBSCRIPT( tmp_subscr_target_4, tmp_subscr_subscript_4 );
    if ( tmp_binop_right_5 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_1 );
        Py_DECREF( tmp_binop_left_5 );

        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }
    tmp_binop_left_4 = BINARY_OPERATION_SUB( tmp_binop_left_5, tmp_binop_right_5 );
    Py_DECREF( tmp_binop_left_5 );
    Py_DECREF( tmp_binop_right_5 );
    if ( tmp_binop_left_4 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_1 );

        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }
    tmp_binop_right_4 = const_int_pos_2;
    tmp_binop_right_1 = POWER_OPERATION( tmp_binop_left_4, tmp_binop_right_4 );
    Py_DECREF( tmp_binop_left_4 );
    if ( tmp_binop_right_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_binop_left_1 );

        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }
    tmp_call_arg_element_1 = BINARY_OPERATION_ADD( tmp_binop_left_1, tmp_binop_right_1 );
    Py_DECREF( tmp_binop_left_1 );
    Py_DECREF( tmp_binop_right_1 );
    if ( tmp_call_arg_element_1 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }
    frame_function->f_lineno = 151;
    tmp_return_value = CALL_FUNCTION_WITH_ARGS1( tmp_called_1, tmp_call_arg_element_1 );
    Py_DECREF( tmp_call_arg_element_1 );
    if ( tmp_return_value == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_function->f_lineno = 151;
        goto frame_exception_exit_1;
    }
    goto frame_return_exit_1;

#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif
    // Put the previous frame back on top.
    popFrameStack();
#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );
    goto frame_no_exception_1;
    frame_return_exit_1:;
#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif
    popFrameStack();
#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );
    goto function_return_exit;
    frame_exception_exit_1:;
#if 0
    RESTORE_FRAME_EXCEPTION( frame_function );
#endif

    if ( exception_tb == NULL )
    {
        exception_tb = MAKE_TRACEBACK( INCREASE_REFCOUNT( frame_function ) );
    }
    else if ( exception_tb->tb_frame != frame_function )
    {
        PyTracebackObject *traceback_new = (PyTracebackObject *)MAKE_TRACEBACK( INCREASE_REFCOUNT( frame_function ) );
        traceback_new->tb_next = exception_tb;
        exception_tb = traceback_new;
    }


    tmp_frame_locals = PyDict_New();
    if ((par_a.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_a,
            par_a.object
        );

    }
    if ((par_b.object != NULL))
    {
        PyDict_SetItem(
            tmp_frame_locals,
            const_str_plain_b,
            par_b.object
        );

    }
    detachFrame( exception_tb, tmp_frame_locals );


    popFrameStack();

#if PYTHON_VERSION >= 340
    frame_function->f_executing -= 1;
#endif
    Py_DECREF( frame_function );

    // Return the error.
    goto function_exception_exit;
    frame_no_exception_1:;


    // Return statement must be present.
    assert(false);
function_exception_exit:
    assert( exception_type );
    PyErr_Restore( exception_type, exception_value, (PyObject *)exception_tb );
    return NULL;
function_return_exit:
    return tmp_return_value;

}
static PyObject *fparse_function_5_mag_difference_of_module_helpers( Nuitka_FunctionObject *self, PyObject **args, Py_ssize_t args_size, PyObject *kw )
{
    assert( kw == NULL || PyDict_Check( kw ) );

    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_size = kw ? PyDict_Size( kw ) : 0;
    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_found = 0;
    NUITKA_MAY_BE_UNUSED Py_ssize_t kw_only_found = 0;
    Py_ssize_t args_given = args_size;
    PyObject *_python_par_a = NULL;
    PyObject *_python_par_b = NULL;
    // Copy given dictionary values to the the respective variables:
    if ( kw_size > 0 )
    {
        Py_ssize_t ppos = 0;
        PyObject *key, *value;

        while( PyDict_Next( kw, &ppos, &key, &value ) )
        {
#if PYTHON_VERSION < 300
            if (unlikely( !PyString_Check( key ) && !PyUnicode_Check( key ) ))
#else
            if (unlikely( !PyUnicode_Check( key ) ))
#endif
            {
                PyErr_Format( PyExc_TypeError, "mag_difference() keywords must be strings" );
                goto error_exit;
            }

            NUITKA_MAY_BE_UNUSED bool found = false;

            Py_INCREF( key );
            Py_INCREF( value );

            // Quick path, could be our value.
            if ( found == false && const_str_plain_a == key )
            {
                assert( _python_par_a == NULL );
                _python_par_a = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && const_str_plain_b == key )
            {
                assert( _python_par_b == NULL );
                _python_par_b = value;

                found = true;
                kw_found += 1;
            }

            // Slow path, compare against all parameter names.
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_a, key ) == 1 )
            {
                assert( _python_par_a == NULL );
                _python_par_a = value;

                found = true;
                kw_found += 1;
            }
            if ( found == false && RICH_COMPARE_BOOL_EQ( const_str_plain_b, key ) == 1 )
            {
                assert( _python_par_b == NULL );
                _python_par_b = value;

                found = true;
                kw_found += 1;
            }


            Py_DECREF( key );

            if ( found == false )
            {
               Py_DECREF( value );

               PyErr_Format(
                   PyExc_TypeError,
                   "mag_difference() got an unexpected keyword argument '%s'",
                   Nuitka_String_Check( key ) ? Nuitka_String_AsString( key ) : "<non-string>"
               );

               goto error_exit;
            }
        }

#if PYTHON_VERSION < 300
        assert( kw_found == kw_size );
        assert( kw_only_found == 0 );
#endif
    }

    // Check if too many arguments were given in case of non star args
    if (unlikely( args_given > 2 ))
    {
#if PYTHON_VERSION < 270
        ERROR_TOO_MANY_ARGUMENTS( self, args_given, kw_size );
#elif PYTHON_VERSION < 330
        ERROR_TOO_MANY_ARGUMENTS( self, args_given + kw_found );
#else
        ERROR_TOO_MANY_ARGUMENTS( self, args_given, kw_only_found );
#endif
        goto error_exit;
    }


    // Copy normal parameter values given as part of the args list to the respective variables:

    if (likely( 0 < args_given ))
    {
         if (unlikely( _python_par_a != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 0 );
             goto error_exit;
         }

        _python_par_a = INCREASE_REFCOUNT( args[ 0 ] );
    }
    else if ( _python_par_a == NULL )
    {
        if ( 0 + self->m_defaults_given >= 2  )
        {
            _python_par_a = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 0 - 2 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }
    if (likely( 1 < args_given ))
    {
         if (unlikely( _python_par_b != NULL ))
         {
             ERROR_MULTIPLE_VALUES( self, 1 );
             goto error_exit;
         }

        _python_par_b = INCREASE_REFCOUNT( args[ 1 ] );
    }
    else if ( _python_par_b == NULL )
    {
        if ( 1 + self->m_defaults_given >= 2  )
        {
            _python_par_b = INCREASE_REFCOUNT( PyTuple_GET_ITEM( self->m_defaults, self->m_defaults_given + 1 - 2 ) );
        }
#if PYTHON_VERSION < 330
        else
        {
#if PYTHON_VERSION < 270
            ERROR_TOO_FEW_ARGUMENTS( self, kw_size, args_given + kw_found );
#elif PYTHON_VERSION < 300
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found );
#else
            ERROR_TOO_FEW_ARGUMENTS( self, args_given + kw_found - kw_only_found );
#endif

            goto error_exit;
        }
#endif
    }

#if PYTHON_VERSION >= 330
    if (unlikely( _python_par_a == NULL || _python_par_b == NULL ))
    {
        PyObject *values[] = { _python_par_a, _python_par_b };
        ERROR_TOO_FEW_ARGUMENTS( self, values );

        goto error_exit;
    }
#endif


    return impl_function_5_mag_difference_of_module_helpers( self, _python_par_a, _python_par_b );

error_exit:;

    Py_XDECREF( _python_par_a );
    Py_XDECREF( _python_par_b );

    return NULL;
}

static PyObject *dparse_function_5_mag_difference_of_module_helpers( Nuitka_FunctionObject *self, PyObject **args, int size )
{
    if ( size == 2 )
    {
        return impl_function_5_mag_difference_of_module_helpers( self, INCREASE_REFCOUNT( args[ 0 ] ), INCREASE_REFCOUNT( args[ 1 ] ) );
    }
    else
    {
        PyObject *result = fparse_function_5_mag_difference_of_module_helpers( self, args, size, NULL );
        return result;
    }

}




static PyObject *MAKE_FUNCTION_function_1_find_crossings_of_module_helpers(  )
{
    PyObject *result = Nuitka_Function_New(
        fparse_function_1_find_crossings_of_module_helpers,
        dparse_function_1_find_crossings_of_module_helpers,
        const_str_plain_find_crossings,
#if PYTHON_VERSION >= 330
        NULL,
#endif
        codeobj_ca92ac72b7622f4242ad42d6ae97b5d7,
        NULL,
#if PYTHON_VERSION >= 300
        NULL,
        const_dict_empty,
#endif
        module_helpers,
        const_str_digest_4524f30e0c6d9d9f108a2d6de2319892
    );

    return result;
}



static PyObject *MAKE_FUNCTION_function_2_do_vectors_intersect_of_module_helpers(  )
{
    PyObject *result = Nuitka_Function_New(
        fparse_function_2_do_vectors_intersect_of_module_helpers,
        dparse_function_2_do_vectors_intersect_of_module_helpers,
        const_str_plain_do_vectors_intersect,
#if PYTHON_VERSION >= 330
        NULL,
#endif
        codeobj_7df00fdabc49ad9593c404a7f3e69637,
        NULL,
#if PYTHON_VERSION >= 300
        NULL,
        const_dict_empty,
#endif
        module_helpers,
        const_str_digest_3a2dc88a511759902b6efe905e17a85b
    );

    return result;
}



static PyObject *MAKE_FUNCTION_function_3_cross_product_of_module_helpers(  )
{
    PyObject *result = Nuitka_Function_New(
        fparse_function_3_cross_product_of_module_helpers,
        dparse_function_3_cross_product_of_module_helpers,
        const_str_plain_cross_product,
#if PYTHON_VERSION >= 330
        NULL,
#endif
        codeobj_641b42b44ac395c476fb706d09df8c7a,
        NULL,
#if PYTHON_VERSION >= 300
        NULL,
        const_dict_empty,
#endif
        module_helpers,
        const_str_digest_65401906cd2c069cbdb366a8b9d824a4
    );

    return result;
}



static PyObject *MAKE_FUNCTION_function_4_sign_of_module_helpers(  )
{
    PyObject *result = Nuitka_Function_New(
        fparse_function_4_sign_of_module_helpers,
        dparse_function_4_sign_of_module_helpers,
        const_str_plain_sign,
#if PYTHON_VERSION >= 330
        NULL,
#endif
        codeobj_fd2d6789440b495c075e8467bf2253e1,
        NULL,
#if PYTHON_VERSION >= 300
        NULL,
        const_dict_empty,
#endif
        module_helpers,
        Py_None
    );

    return result;
}



static PyObject *MAKE_FUNCTION_function_5_mag_difference_of_module_helpers(  )
{
    PyObject *result = Nuitka_Function_New(
        fparse_function_5_mag_difference_of_module_helpers,
        dparse_function_5_mag_difference_of_module_helpers,
        const_str_plain_mag_difference,
#if PYTHON_VERSION >= 330
        NULL,
#endif
        codeobj_71cfbdc35b90f7d090c4273357868645,
        NULL,
#if PYTHON_VERSION >= 300
        NULL,
        const_dict_empty,
#endif
        module_helpers,
        const_str_digest_f3dba3f4d453eccf489afa145311fe8b
    );

    return result;
}



#if PYTHON_VERSION >= 300
static struct PyModuleDef mdef_helpers =
{
    PyModuleDef_HEAD_INIT,
    "helpers",   /* m_name */
    NULL,                /* m_doc */
    -1,                  /* m_size */
    NULL,                /* m_methods */
    NULL,                /* m_reload */
    NULL,                /* m_traverse */
    NULL,                /* m_clear */
    NULL,                /* m_free */
  };
#endif

#define _MODULE_UNFREEZER 1

#if _MODULE_UNFREEZER

#include "nuitka/unfreezing.hpp"

// Table for lookup to find "frozen" modules or DLLs, i.e. the ones included in
// or along this binary.
MOD_INIT_DECL( helpers );
static struct Nuitka_MetaPathBasedLoaderEntry meta_path_loader_entries[] =
{
    { (char *)"helpers", MOD_INIT_NAME( helpers ), NUITKA_COMPILED_MODULE },
    { NULL, NULL, 0 }
};

#endif

// The exported interface to CPython. On import of the module, this function
// gets called. It has to have an exact function name, in cases it's a shared
// library export. This is hidden behind the MOD_INIT_DECL.

MOD_INIT_DECL( helpers )
{

#if defined(_NUITKA_EXE) || PYTHON_VERSION >= 300
    static bool _init_done = false;

    // Packages can be imported recursively in deep executables.
    if ( _init_done )
    {
        return MOD_RETURN_VALUE( module_helpers );
    }
    else
    {
        _init_done = true;
    }
#endif

#ifdef _NUITKA_MODULE
    // In case of a stand alone extension module, need to call initialization
    // the init here because that's the first and only time we are going to get
    // called here.

    // Initialize the constant values used.
    _initBuiltinModule();
    _initConstants();

    // Initialize the compiled types of Nuitka.
    PyType_Ready( &Nuitka_Generator_Type );
    PyType_Ready( &Nuitka_Function_Type );
    PyType_Ready( &Nuitka_Method_Type );
    PyType_Ready( &Nuitka_Frame_Type );
#if PYTHON_VERSION < 300
    initSlotCompare();
#endif

    patchBuiltinModule();
    patchTypeComparison();

#endif

#if _MODULE_UNFREEZER
    registerMetaPathBasedUnfreezer( meta_path_loader_entries );
#endif

    _initModuleConstants();
    _initModuleCodeObjects();

    // puts( "in inithelpers" );

    // Create the module object first. There are no methods initially, all are
    // added dynamically in actual code only.  Also no "__doc__" is initially
    // set at this time, as it could not contain NUL characters this way, they
    // are instead set in early module code.  No "self" for modules, we have no
    // use for it.
#if PYTHON_VERSION < 300
    module_helpers = Py_InitModule4(
        "helpers",       // Module Name
        NULL,                    // No methods initially, all are added
                                 // dynamically in actual module code only.
        NULL,                    // No __doc__ is initially set, as it could
                                 // not contain NUL this way, added early in
                                 // actual code.
        NULL,                    // No self for modules, we don't use it.
        PYTHON_API_VERSION
    );
#else
    module_helpers = PyModule_Create( &mdef_helpers );
#endif

    moduledict_helpers = (PyDictObject *)((PyModuleObject *)module_helpers)->md_dict;

    assertObject( module_helpers );

// Seems to work for Python2.7 out of the box, but for Python3, the module
// doesn't automatically enter "sys.modules", so do it manually.
#if PYTHON_VERSION >= 300
    {
        int r = PyObject_SetItem( PySys_GetObject( (char *)"modules" ), const_str_plain_helpers, module_helpers );

        assert( r != -1 );
    }
#endif

    // For deep importing of a module we need to have "__builtins__", so we set
    // it ourselves in the same way than CPython does. Note: This must be done
    // before the frame object is allocated, or else it may fail.

    PyObject *module_dict = PyModule_GetDict( module_helpers );

    if ( PyDict_GetItem( module_dict, const_str_plain___builtins__ ) == NULL )
    {
        PyObject *value = (PyObject *)builtin_module;

        // Check if main module, not a dict then.
#if !defined(_NUITKA_EXE) || !0
        value = PyModule_GetDict( value );
#endif

#ifndef __NUITKA_NO_ASSERT__
        int res =
#endif
            PyDict_SetItem( module_dict, const_str_plain___builtins__, value );

        assert( res == 0 );
    }

#if PYTHON_VERSION >= 330
#if _MODULE_UNFREEZER
    PyDict_SetItem( module_dict, const_str_plain___loader__, metapath_based_loader );
#else
    PyDict_SetItem( module_dict, const_str_plain___loader__, Py_None );
#endif
#endif

    // Temp variables if any
    PyObject *exception_type, *exception_value;
    PyTracebackObject *exception_tb;
    PyObject *tmp_assign_source_1;
    PyObject *tmp_assign_source_2;
    PyObject *tmp_assign_source_3;
    PyObject *tmp_assign_source_4;
    PyObject *tmp_assign_source_5;
    PyObject *tmp_assign_source_6;
    PyObject *tmp_assign_source_7;
    PyObject *tmp_assign_source_8;
    PyObject *tmp_assign_source_9;
    PyObject *tmp_assign_source_10;
    PyObject *tmp_import_globals_1;
    PyObject *tmp_source_name_1;
    PyObject *tmp_source_name_2;

    // Module code.
    tmp_assign_source_1 = const_str_digest_b41289c3bcec03f4b67da6ac9f622993;
    UPDATE_STRING_DICT0( moduledict_helpers, (Nuitka_StringObject *)const_str_plain___doc__, tmp_assign_source_1 );
    tmp_assign_source_2 = const_str_digest_e8b4e4b3fe31ecb310658b8a0d0f0c8e;
    UPDATE_STRING_DICT0( moduledict_helpers, (Nuitka_StringObject *)const_str_plain___file__, tmp_assign_source_2 );
    // Frame without reuse.
    PyFrameObject *frame_module = MAKE_FRAME( codeobj_c02aa247315f8c1fd9f538dca91c1671, module_helpers );

    // Push the new frame as the currently active one, and we should be exlusively
    // owning it.
    pushFrameStack( frame_module );
    assert( Py_REFCNT( frame_module ) == 1 );

#if PYTHON_VERSION >= 340
    frame_module->f_executing += 1;
#endif

    // Framed code:
    tmp_import_globals_1 = ((PyModuleObject *)module_helpers)->md_dict;
    frame_module->f_lineno = 5;
    tmp_assign_source_3 = IMPORT_MODULE( const_str_plain_numpy, tmp_import_globals_1, tmp_import_globals_1, Py_None, const_int_neg_1 );
    if ( tmp_assign_source_3 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_module->f_lineno = 5;
        goto frame_exception_exit_1;
    }
    UPDATE_STRING_DICT1( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_n, tmp_assign_source_3 );
    tmp_source_name_1 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_n );

    if (unlikely( tmp_source_name_1 == NULL ))
    {
        tmp_source_name_1 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_n );
    }

    if ( tmp_source_name_1 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 0 ], 23, 0 );
        exception_tb = NULL;

        frame_module->f_lineno = 7;
        goto frame_exception_exit_1;
    }

    tmp_assign_source_4 = LOOKUP_ATTRIBUTE( tmp_source_name_1, const_str_plain_sqrt );
    if ( tmp_assign_source_4 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_module->f_lineno = 7;
        goto frame_exception_exit_1;
    }
    UPDATE_STRING_DICT1( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_csqrt, tmp_assign_source_4 );
    tmp_source_name_2 = GET_STRING_DICT_VALUE( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_n );

    if (unlikely( tmp_source_name_2 == NULL ))
    {
        tmp_source_name_2 = GET_STRING_DICT_VALUE( dict_builtin, (Nuitka_StringObject *)const_str_plain_n );
    }

    if ( tmp_source_name_2 == NULL )
    {

        exception_type = INCREASE_REFCOUNT( PyExc_NameError );
        exception_value = UNSTREAM_STRING( &constant_bin[ 0 ], 23, 0 );
        exception_tb = NULL;

        frame_module->f_lineno = 8;
        goto frame_exception_exit_1;
    }

    tmp_assign_source_5 = LOOKUP_ATTRIBUTE( tmp_source_name_2, const_str_plain_floor );
    if ( tmp_assign_source_5 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );


        frame_module->f_lineno = 8;
        goto frame_exception_exit_1;
    }
    UPDATE_STRING_DICT1( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_floor, tmp_assign_source_5 );
    tmp_assign_source_6 = MAKE_FUNCTION_function_1_find_crossings_of_module_helpers(  );
    if ( tmp_assign_source_6 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_assign_source_6 );

        frame_module->f_lineno = 11;
        goto frame_exception_exit_1;
    }
    UPDATE_STRING_DICT1( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_find_crossings, tmp_assign_source_6 );
    tmp_assign_source_7 = MAKE_FUNCTION_function_2_do_vectors_intersect_of_module_helpers(  );
    if ( tmp_assign_source_7 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_assign_source_7 );

        frame_module->f_lineno = 121;
        goto frame_exception_exit_1;
    }
    UPDATE_STRING_DICT1( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_do_vectors_intersect, tmp_assign_source_7 );
    tmp_assign_source_8 = MAKE_FUNCTION_function_3_cross_product_of_module_helpers(  );
    if ( tmp_assign_source_8 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_assign_source_8 );

        frame_module->f_lineno = 140;
        goto frame_exception_exit_1;
    }
    UPDATE_STRING_DICT1( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_cross_product, tmp_assign_source_8 );
    tmp_assign_source_9 = MAKE_FUNCTION_function_4_sign_of_module_helpers(  );
    if ( tmp_assign_source_9 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_assign_source_9 );

        frame_module->f_lineno = 145;
        goto frame_exception_exit_1;
    }
    UPDATE_STRING_DICT1( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_sign, tmp_assign_source_9 );
    tmp_assign_source_10 = MAKE_FUNCTION_function_5_mag_difference_of_module_helpers(  );
    if ( tmp_assign_source_10 == NULL )
    {
        assert( ERROR_OCCURED() );

        PyErr_Fetch( &exception_type, &exception_value, (PyObject **)&exception_tb );
        Py_DECREF( tmp_assign_source_10 );

        frame_module->f_lineno = 149;
        goto frame_exception_exit_1;
    }
    UPDATE_STRING_DICT1( moduledict_helpers, (Nuitka_StringObject *)const_str_plain_mag_difference, tmp_assign_source_10 );

    // Restore frame exception if necessary.
#if 0
    RESTORE_FRAME_EXCEPTION( frame_module );
#endif
    popFrameStack();

    assertFrameObject( frame_module );
    Py_DECREF( frame_module );

    goto frame_no_exception_1;
    frame_exception_exit_1:;
#if 0
    RESTORE_FRAME_EXCEPTION( frame_module );
#endif

    if ( exception_tb == NULL )
    {
        exception_tb = MAKE_TRACEBACK( INCREASE_REFCOUNT( frame_module ) );
    }
    else if ( exception_tb->tb_frame != frame_module )
    {
        PyTracebackObject *traceback_new = (PyTracebackObject *)MAKE_TRACEBACK( INCREASE_REFCOUNT( frame_module ) );
        traceback_new->tb_next = exception_tb;
        exception_tb = traceback_new;
    }

    // Put the previous frame back on top.
    popFrameStack();

#if PYTHON_VERSION >= 340
    frame_module->f_executing -= 1;
#endif
    Py_DECREF( frame_module );

    // Return the error.
    goto module_exception_exit;
    frame_no_exception_1:;

    return MOD_RETURN_VALUE( module_helpers );
module_exception_exit:
    PyErr_Restore( exception_type, exception_value, (PyObject *)exception_tb );
    return MOD_RETURN_VALUE( NULL );
}
