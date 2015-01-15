#include <iostream>
#include <map>
#include "H_BeamParticle.h"
using namespace std;

void filterFP420(){
  unsigned int line;
  H_BeamParticle * part;
  std::map< unsigned int, H_BeamParticle* >::iterator it;
  
  bool is_stop;
  int direction;
  
  float x1_420;
  float y1_420;
  
  if ( m_beamPart.size() && lengthfp420>0. ) {
    
    for (it = m_beamPart.begin(); it != m_beamPart.end(); ++it ) {
      line = (*it).first;
      part = (*it).second;
      
      if ( (*m_isCharged.find( line )).second ) {
	direction = (*m_direct.find( line )).second;
	if ( direction == 1 ) {
	  part->computePath( m_beamlineFP4201, 1 );
	  is_stop = part->stopped( m_beamlineFP4201 );
	}
	else if ( direction == -1 ){
	  part->computePath( m_beamlineFP4202, -1 );
	  is_stop = part->stopped( m_beamlineFP4202 );
	}
	
	//propagating
	m_isStoppedfp420[line] = is_stop;
	
	if (!is_stop) {
	  if ( direction == 1 ) part->propagate( m_rpp420_f );
	  if ( direction == -1 ) part->propagate( m_rpp420_b );
	  x1_420 = part->getX();
	  y1_420 = part->getY();
	  if(m_verbosity) cout << "=== Hector:filterFP420: x=  "<< x1_420 <<" y= " << y1_420 << endl;
	  
	  m_xAtTrPoint[line]  = x1_420;
	  m_yAtTrPoint[line]  = y1_420;
	  m_TxAtTrPoint[line] = part->getTX();
	  m_TyAtTrPoint[line] = part->getTY();
	  m_eAtTrPoint[line]  = part->getE();
	  
	}
      }// if isCharged
      else {
	m_isStoppedfp420[line] = true;// imply that neutral particles stopped to reach 420m
	if(m_verbosity) cout << "=== Hector:filterFP420:  isStopped=" << (*m_isStoppedfp420.find(line)).second <<  endl;
      }
      
    } // for (it = m_beamPart.begin(); it != m_beamPart.end(); it++ ) 
  } // if ( m_beamPart.size() )
  
}

int main() {

 return 0;
}

