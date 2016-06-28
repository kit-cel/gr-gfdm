/* -*- c++ -*- */

#define GFDM_API
#define DIGITAL_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "gfdm_swig_doc.i"

%{
#include "gfdm/framer_cc.h"
#include "gfdm/modulator_cc.h"
#include "gfdm/receiver_cc.h"
#include "gfdm/advanced_receiver_cc.h"
#include "gnuradio/digital/constellation.h"
#include "gfdm/sync_cc.h"
#include "gfdm/cyclic_prefixer_cc.h"
#include "gfdm/preamble_generator.h"
#include "gfdm/remove_prefix_cc.h"
#include "gfdm/simple_modulator_cc.h"
#include "gfdm/modulator_kernel_cc.h"
#include "gfdm/add_cyclic_prefix_cc.h"
#include "gfdm/simple_receiver_cc.h"
%}

%include "gfdm/framer_cc.h"
%include "gfdm/modulator_cc.h"
%include "gfdm/receiver_cc.h"
%include "gfdm/advanced_receiver_cc.h"
%include "gnuradio/digital/constellation.h"

GR_SWIG_BLOCK_MAGIC2(gfdm, framer_cc);
GR_SWIG_BLOCK_MAGIC2(gfdm, modulator_cc);
GR_SWIG_BLOCK_MAGIC2(gfdm, receiver_cc);
GR_SWIG_BLOCK_MAGIC2(gfdm, advanced_receiver_cc);

%include "gnuradio/swig/constellation.i"
%include "gfdm/sync_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, sync_cc);
%include "gfdm/cyclic_prefixer_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, cyclic_prefixer_cc);
%include "gfdm/preamble_generator.h"
%include "preamble_generator.i"
%include "gfdm/remove_prefix_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, remove_prefix_cc);
%include "gfdm/simple_modulator_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, simple_modulator_cc);
%include "gfdm/simple_receiver_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, simple_receiver_cc);
