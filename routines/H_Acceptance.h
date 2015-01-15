#ifndef _H_Acceptance_
#define _H_Acceptance_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_Acceptance.h
/// \brief Collection of routines about physics output from Hector, excepting explicit variable reconstruction (see also H_Reconstruction.cpp)

// C++ includes
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>

// local includes
#include "H_AbstractBeamLine.h"
#include "H_Beam.h"
#include "H_BeamLine.h"
#include "H_BeamParticle.h"
#include "H_OpticalElement.h"
#include "H_RomanPot.h"


// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TEllipse.h"
using namespace std;

/// Draws the 1D-acceptance for any optical element, for a given energy range and a given virtuality
TH1F* acceptance_element(H_AbstractBeamLine * beamline, const string element_name, const unsigned int number_of_particles=10, const float E_min=0., const float E_max=150., const unsigned int Ne = 10, const float Q2 =-0.1, const bool NonLinear=true, const int crang_sign=-1);
	/// @param beamline : beamline to test
	/// @param number_of_particles : in the beam
	/// @param E_min : min value for the beam particle energy loss
	/// @param E_max : max value for  the beam particle energy loss
	/// @param Ne : number of steps in energy
	/// @param Q2 : virtuality value (t) of the emitted gamma
	/// @param NonLinear : (des)activates the NonLinear effects for the beam propagation. Should be always on.


/// Draws the lego plot of irradiation levels, from \f$ pp \rightarrow pX \f$
void acceptance_fluence(float rpos = 220., float deltax = 2000., string file = "data/LHCB1IR5_v6.500.tfs", int side = 1, const int N = 1000, const bool save=false);
	/// @param rpos is the s coordinate where to look at the irradiation, in [m] from IP
	/// @param deltax is the x coordinate where a fake roman pot is set, in [microm] from the beam center
	/// @param file is the optics table
	/// @param side is the direction of propagation (1 is downstream, -1 is upstream)
	/// @param N is the number of particles to propagate (ideal one : 10^6)
	/// @param save : boolean

/// Returns the usual beamline with a RP
H_AbstractBeamLine * createBeamline(float rp_pos, float rp_x, string filename, int side, string& rp_name);
	/// @param rp_pos is the RP position in [m], from the IP
	/// @param rp_x is the RP position in [\f$ \mu m \f$], when the nominal beam is assumed to be in 0.
	/// @param filename : optics source file
	/// @param side is the direction of propagation (1 : forward, -1 : backward)
	/// @param rp_name : returns the RP name, for subsequent use of this specific optical element


/// Draws the acceptance for each element of the beamline hit by some particles, from the IP until the RP.
void acceptance_beamline(float rp_pos, float rp_x, string filename="data/LHCB1IR5_v6.500.tfs", int side=1, int crang_sign=-1, float phimin = 0, float phimax = 2*3.141592653589, const float log10x_min=-3.5, const float log10x_max=-0.3, const unsigned int Ne = 20, const float log10t_min=-3, const float log10t_max=0, const unsigned int Nq = 33, unsigned int number_of_particles=30, bool NonLinear = true, bool draw_RP_only=false); 
/// @param rp_pos : distance of the roman pots, in meters from the IP
/// @param rp_x   : roman pot transverse position from the center of the beam, in micrometers [\f$ \mu m \f$]
/// @param filename : name of the optics source file
/// @param side : direction of propagation (+1 : forward, -1: backward)
/// @param crang_sign : sign of the (half) crossign angle (should be -1 when forward)
/// @param phimin : lower bound for the random phi angle, in H_Particle::emitGamma()
/// @param phimax : higher bound for the random phi angle, in H_Particle::emitGamma()
/// @param Ne : number of steps in energy
/// @param log10x_min : min value for x (\f$ x = \xi \f$)
/// @param log10x_max : max value for x (\f$ x = \xi \f$)
/// @param Nq : number of steps in t
/// @param log10t_min : min of the virtuality axis / t value (\f$ t = -Q^2 \f$)
/// @param log10t_max : max of the virtuality axis / t value (\f$ t = -Q^2 \f$)
/// @param number_of_particles : in the beam
/// @param NonLinear : (des)activates the NonLinear effects for the beam propagation. Should be always on.
/// @param draw_RP_only : draw option


/// Same routine as acceptance_beamline, but dedicated to RP only. Should be slightly faster.
void acceptance_rp(float rp_pos, float rp_x, string filename="data/LHCB1IR5_v6.500.tfs", const int side=1, const int crang_sign = -1, unsigned int number_of_particles=10, const float E_min=0, const float E_max=3000, const unsigned int Ne = 10, const float log10t_min=-3, const float log10t_max=0, const unsigned int Nq = 10, bool NonLinear=true, const char * add2title="", string outfilename="");
/// @param rp_pos : distance of the roman pots, in meters from the IP
/// @param rp_x   : roman pot transverse position from the center of the beam, in micrometers [\f$\mu m\f$]
/// @param filename : name of the optics source file
/// @param side : direction of propagation (+1 : forward, -1: backward)
/// @param crang_sign : sign of the (half) crossign angle (should be -1 when forward)
/// @param number_of_particles : in the beam
/// @param E_min : min value for the beam particle energy loss
/// @param E_max : max value for  the beam particle energy loss
/// @param Ne : number of steps in energy
/// @param log10t_min : min of the virtuality axis / t value (\f$ t = -Q^2 \f$)
/// @param log10t_max : max of the virtuality axis / t value (\f$ t = -Q^2 \f$)
/// @param Nq : number of steps in t
/// @param NonLinear : (des)activates the NonLinear effects for the beam propagation. Should be always on.
/// @param add2title : extra string added to the histogram title
/// @param outfilename : name of file to save to


/// Draws the 1D-acceptance for a RP detector, for two given virtuality values
void acceptance_rp_1D(float rp_pos=420., float rp_x=4000., string filename="data/LHCB1IR5_v6.500.tfs", const int side=1, const int crang_sign = -1, unsigned int number_of_particles=10, const float E_min=0., const float E_max=150., const unsigned int Ne = 10, const float Q2_0 =-0.1, const float Q2_1=-1, bool NonLinear=true, const char * add2title="", string outfilename=""); 
/// @param rp_pos : distance of the roman pots, in meters from the IP
/// @param rp_x   : roman pot transverse position from the center of the beam, in micrometers [\f$\mu m\f$]
/// @param filename : name of the optics source file
/// @param side : direction of propagation (+1 : forward, -1: backward)
/// @param crang_sign : sign of the (half) crossign angle (should be -1 when forward)
/// @param number_of_particles : in the beam
/// @param E_min : min value for the beam particle energy loss
/// @param E_max : max value for  the beam particle energy loss
/// @param Ne : number of steps in energy
/// @param Q2_0 : first virtuality value
/// @param Q2_1 : second virtuality value
/// @param NonLinear : (des)activates the NonLinear effects for the beam propagation. Should be always on.
/// @param add2title : extra string added to the histogram title
/// @param outfilename : name of file to save to


/// draws the aperture shape of a given optical element, as well as protons "sensing" this aperture
void acceptance_profile(string element_name = "\"MB.B9R5.B1\"") ;
	/// @param element_name is the id of the optical element whose aperture is being drawn

#endif

