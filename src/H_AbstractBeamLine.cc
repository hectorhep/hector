  /* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                         *
*                   --<--<--  A fast simulator --<--<--     *
*                 / --<--<--     of particle   --<--<--     *
*  ----HECTOR----<                                          *
*                 \ -->-->-- transport through -->-->--     *
*                   -->-->-- generic beamlines -->-->--     *
*                                                           *
* JINST 2:P09005 (2007)                                     *
*      X Rouby, J de Favereau, K Piotrzkowski (CP3)         *
*       http://www.fynu.ucl.ac.be/hector.html               *
*                                                           *
* Center for Cosmology, Particle Physics and Phenomenology  *
*              Universite catholique de Louvain             *
*                 Louvain-la-Neuve, Belgium                 *
 *                                                         *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/// \file H_AbstractBeamLine.cc
/// \brief Class describing ideal beamline.
///
/// Units : angles [rad], distances [m], energies [GeV], c=[1].

// c++ #includes
#include <iostream>
#include <cmath>

// local #includes
#include "H_Parameters.h"
#include "H_TransportMatrices.h"
#include "H_Drift.h"
#include "H_AbstractBeamLine.h"
#include "H_BeamParticle.h"
#include "H_RomanPot.h"
using namespace std;

void H_AbstractBeamLine::init(const float length, const float energy) {
	beam_mat.ResizeTo(MDIM,MDIM);
	beam_mat = driftmat(length);
  	beam_length = length;
	beam_energy = energy;
	H_Drift * drift0 = new H_Drift("Drift0",0.,length,energy);
	add(drift0);
	return;
}

H_AbstractBeamLine::H_AbstractBeamLine(const H_AbstractBeamLine& beamline) : 
	matrices(beamline.matrices)   {
	//elements = beamline.elements; //<-- bad ! the new vector contains the same pointers as the previous one
	cloneElements(beamline);
	beam_mat.ResizeTo(MDIM,MDIM);
	beam_mat = beamline.beam_mat;
	beam_length = beamline.beam_length;
	beam_energy = beamline.beam_energy;
}

H_AbstractBeamLine& H_AbstractBeamLine::operator=(const H_AbstractBeamLine& beamline) {
	if(this== &beamline) return *this;
	//elements = beamline.elements; //<-- bad ! the new vector contains the same pointers as the previous one
	cloneElements(beamline);
        matrices = beamline.matrices;
        beam_mat = beamline.beam_mat;
	beam_length = beamline.beam_length;
	beam_energy = beamline.beam_energy;
	return *this;
}

H_AbstractBeamLine* H_AbstractBeamLine::clone() const {
	H_AbstractBeamLine* temp_beam = new H_AbstractBeamLine(beam_length, beam_energy);
	vector<H_OpticalElement*>::const_iterator element_i;
	for (element_i = elements.begin(); element_i<elements.end(); element_i++) {
		if((*element_i)->getType()!=DRIFT) {
			H_OpticalElement* temp_el = (*element_i)->clone();
			temp_beam->add(temp_el);
		}
	}
	temp_beam->beam_mat = beam_mat;
	temp_beam->matrices = matrices;
	return temp_beam;
}

H_AbstractBeamLine::~H_AbstractBeamLine() {
	vector<H_OpticalElement*>::iterator element_i;
	for (element_i = elements.begin(); element_i<elements.end(); element_i++) {
		delete (*element_i);
	}
	elements.clear(); 
	matrices.clear(); 
}

void H_AbstractBeamLine::cloneElements(const H_AbstractBeamLine& beam) {
	vector<H_OpticalElement*>::const_iterator element_i;
        for (element_i = beam.elements.begin(); element_i< beam.elements.end(); element_i++) {
        	H_OpticalElement* temp_el = (*element_i)->clone();
                elements.push_back(temp_el);
        }
}

void H_AbstractBeamLine::add(H_OpticalElement * newElement) {
	/// @param newElement is added to the beamline
//	H_OpticalElement * el = new H_OpticalElement(*newElement);
//	H_OpticalElement * el = const_cast<H_OpticalElement*> newElement; 
//	H_OpticalElement * el = newElement;
	elements.push_back(newElement);
 	float a = newElement->getS()+newElement->getLength();
	if (a > beam_length)	{
		beam_length = a;
		if(VERBOSE) cout<<"<H_AbstractBeamLine> WARNING : element ("<< newElement->getName()<<") too far away. The beam length has been extended to "<< beam_length << ". "<<endl;
	}
	calcSequence();
	calcMatrix();
	return;
}

const TMatrix H_AbstractBeamLine::getBeamMatrix() const {
	return beam_mat;
}

const TMatrix H_AbstractBeamLine::getBeamMatrix(const float eloss,const float p_mass, const float p_charge) {

	vector<H_OpticalElement*>::iterator element_i;
    TMatrix calc_mat(MDIM,MDIM);

    // initialization
	calc_mat.UnitMatrix();
	
	// multiplies the matrix of each beam's element
	// and add each product matrix to the list of matrices.
	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
		calc_mat *= (*element_i)->getMatrix(eloss,p_mass,p_charge);
	}
	return calc_mat;
}

const TMatrix H_AbstractBeamLine::getPartialMatrix(const string& elname, const float eloss, const float p_mass, const float p_charge) {

	vector<H_OpticalElement*>::iterator element_i;
	TMatrix calc_mat(MDIM,MDIM);

	calc_mat.UnitMatrix();

	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
		calc_mat *= (*element_i)->getMatrix(eloss,p_mass,p_charge);
		if(elname==(*element_i)->getName()) {
			return calc_mat;
		}
	}
	cout<<"<H_AbstractBeamLine> Element "<<elname<<" desn't exist. Returning full beam matrix"<<endl;
	return calc_mat;
}

const TMatrix H_AbstractBeamLine::getPartialMatrix(const unsigned int element_position) const {
 	//const int N = (element_position<0)?0:(( (element_position)>elements.size()-1)?elements.size()-1:element_position);
	const int N = (element_position>elements.size()-1)?elements.size()-1:element_position;
	return *(matrices.begin()+N); // //for optimization of the code :same as return &matrices[N];
}

const TMatrix H_AbstractBeamLine::getPartialMatrix(const H_OpticalElement * element) const{
	// returns the transport matrix to transport until the end of the specified element
	// !!! 2 elements should never have the same name in "elements" !!!

	vector<H_OpticalElement*>::const_iterator element_i;
	vector<TMatrix>::const_iterator matrix_i;
	TMatrix calc_mat(MDIM,MDIM);

	// parses the list of optical elements and find the searched one
	for(element_i = elements.begin(),matrix_i = matrices.begin(); element_i < elements.end(); element_i++, matrix_i++) {
		if(element->getName() == (*element_i)->getName()) {
			// element has been found
			calc_mat = *matrix_i;
        	}
	}
	return calc_mat;
}

H_OpticalElement * H_AbstractBeamLine::getElement(const unsigned int element_position) {
	const unsigned int N = (element_position>elements.size())?elements.size():element_position;
	return *(elements.begin()+N);//for optimization of the code :same as return &elements[N];
}

H_OpticalElement * H_AbstractBeamLine::getElement(const unsigned int element_position) const {
	const unsigned int N = (element_position>elements.size())?elements.size():element_position;
        return *(elements.begin()+N);//for optimization of the code :same as return &elements[N];
}


H_OpticalElement * H_AbstractBeamLine::getElement(const string& el_name) {
	for(unsigned int i=0; i < elements.size(); i++) {
		if( (*(elements.begin()+i))->getName() == el_name ) 
			return *(elements.begin()+i);
	} // if found -> return ; else : not found at all !	
	cout<<"<H_AbstractBeamLine> Element "<<el_name<<" not found"<<endl;
	return *(elements.begin()+1);
}

H_OpticalElement * H_AbstractBeamLine::getElement(const string& el_name) const {
        for(unsigned int i=0; i < elements.size(); i++) {
		if( (*(elements.begin()+i))->getName() == el_name)
			return *(elements.begin()+i);
	} // if found -> return ; else : not found at all !
	cout<<"<H_AbstractBeamLine> Element "<<el_name<<" not found"<<endl;
	return *(elements.begin()+1);
}

std::ostream& operator<< (std::ostream& os, const H_AbstractBeamLine& be) {
        vector<H_OpticalElement*>::const_iterator element_i;
        os << "Beamline content" << endl;
        for (element_i = be.elements.begin(); element_i < be.elements.end(); element_i++) {
                os << (int)(element_i - be.elements.begin()) << "\t" << (*element_i)->getName() << "\t" << (*element_i)->getS() << endl;
        }
  return os;
}

void H_AbstractBeamLine::printProperties() const { cout << *this; return; }

void H_AbstractBeamLine::showElements() const{
	vector<H_OpticalElement*>::const_iterator element_i;
	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
		(*element_i)->printProperties();
	}
	cout << "Beam length = " << beam_length << endl;
	cout << "Number of elements (including drifts) = " << getNumberOfElements() << endl;
	return;
}

void H_AbstractBeamLine::showElements(const int type_el) const{
	vector<H_OpticalElement*>::const_iterator element_i;
	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
		if ((*element_i)->getType()==type_el)
			(*element_i)->printProperties();
	}
	return;
}

void H_AbstractBeamLine::showMatrix() const {
	cout << "Transport matrix for the whole beam : " << endl;
	cout << "(x,x',...) = (x*,x'*,...) M " <<endl;
  	printMatrix(beam_mat);
	return;
}

void H_AbstractBeamLine::showMatrices() const{
	// prints the list of all transport matrices, from the whole beam.

	vector<TMatrix>::const_iterator matrix_i;
	vector<H_OpticalElement*>::const_iterator element_i;
	TMatrix temp(MDIM,MDIM);

	for(matrix_i = matrices.begin(), element_i = elements.begin(); matrix_i < matrices.end(); matrix_i++, element_i++) {
		temp = *matrix_i;
		cout << "Matrix for transport until s=" << (*element_i)->getS() + (*element_i)->getLength() << "m (" << (*element_i)->getName() << "). " << endl;
		printMatrix(temp);
		cout << endl;
	}
	return ;
}

void H_AbstractBeamLine::calcSequence() {
		// reorders the elements, computes the drifts;

	vector<H_OpticalElement*> temp_elements;
	vector<H_OpticalElement*>::iterator element_i;
		// element_i is a pointer to elements[i]

	if(elements.size()==1) { return; }

	// getting rid of drifts before calculating
	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
		if((*element_i)->getType() == DRIFT) {delete (*element_i); elements.erase(element_i); }
	}

	// ordering the elements in position
	sort(elements.begin(),elements.end(),ordering());
	// inserting the drifts before the other elements
	float current_pos = 0;
	float drift_length=0;

	for(element_i=elements.begin(); element_i < elements.end(); element_i++) {
		drift_length = (*element_i)->getS() - current_pos;
		if(drift_length>0) {
  			H_Drift *dr = new H_Drift(current_pos,drift_length,beam_energy);
			temp_elements.push_back(dr); 
		}
		temp_elements.push_back(*element_i);
		current_pos = (*element_i)->getS() + (*element_i)->getLength();
	}
	
	//adding the last drift
	drift_length = beam_length - current_pos;
	if (drift_length>0) {
			H_Drift *dr = new H_Drift(current_pos,drift_length,beam_energy);
			temp_elements.push_back(dr);
	}
	
	// cleaning : avoid some memory leaks
	elements.clear();

	for(element_i=temp_elements.begin(); element_i < temp_elements.end(); element_i++) {
		elements.push_back(*element_i);
	}
}

void H_AbstractBeamLine::calcMatrix() {
	// computes the transport matrix for the beam upto here...
	vector<H_OpticalElement*>::iterator element_i;
	TMatrix calc_mat(MDIM,MDIM);

	// initialization
	matrices.clear();
	calc_mat.UnitMatrix();

	// multiplies the matrix of each beam's element
	// and add each product matrix to the list of matrices.
	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
		calc_mat *= (*element_i)->getMatrix();
		matrices.push_back(calc_mat);
  	}

	beam_mat.ResizeTo(MDIM,MDIM);
	beam_mat = calc_mat;
	return;
}

float qh(float k) {
        float beta = (log((float)10.0))/0.05;
		// put (std::log((float)10.0)) instead of log(10) to avoid compilation errors
        return 0.8*(1-exp(-beta*fabs(k)));
}

float dh(float k) {
        float psi = (log((float)10.0))/0.002;
		// put (std::log((float)10.0)) instead of log(10) to avoid compilation errors
        return 0.8*(1-exp(-psi*fabs(k)));
}

void H_AbstractBeamLine::draw(const float xmin, const float ymin, const float xmax, const float ymax) const
	return;
}

void H_AbstractBeamLine::drawX(const float a_min, const float a_max, const float scale) const{
	/// @param a_min defines the size of the drawing
	/// @param a_max defines the size of the drawing
	/// @param scale allows to multiply the drawing, i.e. changing the units
	const int N = getNumberOfElements();
	for(int i=0;i<N;i++) {
		float height = fabs(a_max);
        	float meight = fabs(a_min);
        	float size = (height>meight)?meight:height;
		float middle = getElement(i)->getX()*URAD*scale;
        	if(getElement(i)->getType()!=DRIFT) getElement(i)->draw(middle+size/2.,middle-size/2.);
    	}
}

void H_AbstractBeamLine::drawY(const float a_min, const float a_max) const{
        /// @param a_min defines the size of the drawing 
        /// @param a_max defines the size of the drawing
        const int N = getNumberOfElements();
        for(int i=0;i<N;i++) {
                float height = fabs(a_max);
                float meight = fabs(a_min);
                float size = (height>meight)?meight:height;
                float middle = getElement(i)->getY()*URAD;
                if(getElement(i)->getType()!=DRIFT) getElement(i)->draw(middle+size/2.,middle-size/2.);
        }
}

void H_AbstractBeamLine::moveElement(const string& name, const float new_s) {
	/// @param name identifies the element to move
	/// @param new_s is where to put it
       vector<H_OpticalElement*>::iterator element_i;
       for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
	       if(name==(*element_i)->getName()) { (*element_i)->setS(new_s); }
       }

       calcSequence();
       calcMatrix();
       return;
}

void H_AbstractBeamLine::alignElement(const string& name, const float disp_x, const float disp_y) {
	/// @param name identifies the element to move
	/// @param disp_x identifies the displacement to add in x [\f$ \mu m \f$]
	/// @param disp_y identifies the displacement to add in y [\f$ \mu m \f$]
		vector<H_OpticalElement*>::iterator element_i;
		for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
			if(name==(*element_i)->getName()) { 
				(*element_i)->setX((*element_i)->getX()+disp_x);
				(*element_i)->setY((*element_i)->getY()+disp_y);
				return ;
			}
		}
		cout<<"<H_AbstractBeamLine> WARNING : Element "<<name<<" not found."<<endl; 
		return;
}

void H_AbstractBeamLine::tiltElement(const string& name, const float ang_x, const float ang_y) {
	/// @param name identifies the element to move
	/// @param ang_x identifies the angle to add in x
	/// @param ang_y identifies the angle to add in y
		vector<H_OpticalElement*>::iterator element_i;
		for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
			if(name==(*element_i)->getName()) {
				(*element_i)->setTX((*element_i)->getTX()+ang_x);
				(*element_i)->setTY((*element_i)->getTY()+ang_y);
				return ;
			}
		}
		cout<<"<H_AbstractBeamLine> WARNING : Element "<<name<<" not found."<<endl; 
		return;
}

void H_AbstractBeamLine::offsetElements(const float start, const float offset) {
	/// @param start After this s [m] coordinate, all elements will be offset.
	/// @param offset In meters

	extern int relative_energy;
	if(!relative_energy) {
		vector<H_OpticalElement*>::iterator element_i;
   	     for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
			if((*element_i)->getS() > start ) {
    	       (*element_i)->setX(offset);
       	    }
		}
	}
}
