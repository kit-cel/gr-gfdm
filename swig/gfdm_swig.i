/* -*- c++ -*- */

#define GFDM_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "gfdm_swig_doc.i"

%{
#include "gfdm/transmitter_cc.h"
#include "gfdm/framer_cc.h"
//#include "gfdm/receiver_cc.h"
%}


%include "gfdm/transmitter_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, transmitter_cc);
//%include "gfdm/receiver_cc.h"
//GR_SWIG_BLOCK_MAGIC2(gfdm, receiver_cc);
%include "gfdm/framer_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, framer_cc);
