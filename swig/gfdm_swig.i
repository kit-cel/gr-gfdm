/* -*- c++ -*- */

#define GFDM_API
#define DIGITAL_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "gfdm_swig_doc.i"

%{
#include "gfdm/framer_cc.h"
#include "gfdm/modulator_cc.h"
#include "gnuradio/digital/constellation.h"
#include "gfdm/sync_cc.h"
#include "gfdm/cyclic_prefixer_cc.h"
#include "gfdm/preamble_generator.h"
#include "gfdm/remove_prefix_cc.h"
#include "gfdm/simple_modulator_cc.h"
#include "gfdm/modulator_kernel_cc.h"
#include "gfdm/receiver_kernel_cc.h"
#include "gfdm/add_cyclic_prefix_cc.h"
#include "gfdm/improved_sync_algorithm_kernel_cc.h"
#include "gfdm/simple_receiver_cc.h"
#include "gfdm/advanced_receiver_sb_cc.h"
#include "gfdm/advanced_receiver_tsb_cc.h"
#include "gfdm/resource_mapper_kernel_cc.h"
#include "gfdm/resource_mapper_cc.h"
#include "gfdm/frame_energy_detector_cc.h"
#include "gfdm/simple_preamble_sync_cc.h"
#include "gfdm/resource_demapper_cc.h"
#include "gfdm/extract_burst_cc.h"
#include "gfdm/channel_estimator_cc.h"
%}

%include "gfdm/framer_cc.h"
%include "gfdm/modulator_cc.h"
%include "gnuradio/digital/constellation.h"

GR_SWIG_BLOCK_MAGIC2(gfdm, framer_cc);
GR_SWIG_BLOCK_MAGIC2(gfdm, modulator_cc);

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
%include "gfdm/improved_sync_algorithm_kernel_cc.h"
%include "gfdm/advanced_receiver_sb_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, advanced_receiver_sb_cc);
%include "gfdm/advanced_receiver_tsb_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, advanced_receiver_tsb_cc);
%include "gfdm/resource_mapper_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, resource_mapper_cc);
%include "gfdm/frame_energy_detector_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, frame_energy_detector_cc);
%include "gfdm/simple_preamble_sync_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, simple_preamble_sync_cc);
%include "gfdm/resource_demapper_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, resource_demapper_cc);
%include "gfdm/extract_burst_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, extract_burst_cc);
%include "gfdm/channel_estimator_cc.h"
GR_SWIG_BLOCK_MAGIC2(gfdm, channel_estimator_cc);
