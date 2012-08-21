/* Copyright 2008 Guillem Rigaill <guillem.rigaill@curie.fr> 

   This file is part of colibri design for a fast segmentation in the mean

   Colibri is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   Colibri is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with Colibri; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef LISTEH
#define LISTEH
#include "polynome2.h"
#include "numlib.h"


#define INLINED
class Liste {
 private:
     double max, min;
	 Polynome2 *poly;
	 Liste *next;
 public:
    /* constructors and destructors */
#ifdef INLINED
	inline
#endif
    Liste()
		: max(0.), min(0.), next(NULL), poly(NULL) {}
	Liste(double max_, double min_)
		: max(max_), min(min_), next(NULL) , poly(NULL){}
	Liste(double max_, double min_, Polynome2 *poly_)
		: max(max_), min(min_), next(NULL), poly(poly_) {}
	Liste(Polynome2 *poly_)
		: max(0.), min(0.), next(NULL), poly(poly_) {}
	~Liste(){
		delete next;
		delete poly;
	}
	/* fonction setter and getter */
#ifdef INLINED
	inline
#endif
    double getMax();
#ifdef INLINED
	inline
#endif
    void setMax(double max_);
#ifdef INLINED
	inline
#endif
	double getMin();
#ifdef INLINED
	inline
#endif
    void setMin(double min_);

#ifdef INLINED
	inline
#endif
	void setPolynome(Polynome2 * poly_);
#ifdef INLINED
	inline
#endif
	Polynome2 *getPolynome();

#ifdef INLINED
	inline
#endif
	Liste * getNext();
#ifdef INLINED
	inline
#endif
	void setNext(Liste * next_);

	/* Useful */
#ifdef INLINED
	inline
#endif
	void setToNull();
#ifdef INLINED
	inline
#endif
	void insert(Liste * maillon_);
#ifdef INLINED
        inline
#endif
        int compte();
#ifdef INLINED
        inline
#endif
	Liste * removeDoublon();
#ifdef INLINED
        inline
#endif
	void checkForDoublon();

	/* show and others */
	void show();
	void showAllNext();
	
	/* */
#ifdef INLINED
        inline
#endif
	void computeRoots(double);
#ifdef INLINED
        inline
#endif
	void add(double, double, double);
#ifdef INLINED
        inline
#endif
	void computeMinOrMax(double*, int*);
	void resetMaillonBorders(Polynome2*);
	void resetAllBorders(Polynome2*);
	
};

#ifdef INLINED
/* Setter and Getter */
/* */
double Liste::getMax()
{
	return(this->max);
}

void Liste::setMax(double max_)
{
	this->max = max_;
}

double Liste::getMin()
{
	return(this->min);
}

void Liste::setMin(double min_)
{
	this->min = min_;
}

void Liste::setPolynome(Polynome2 * poly_)
{
		this->poly=poly_;
}
Polynome2* Liste::getPolynome()
{
		return(this->poly	);
}
/* */
Liste * Liste::getNext()
{
	return(this->next);
}

void Liste::setNext(Liste * next_)
{
	this->next = next_;
}

void Liste::setToNull()
{
	max=NULL;
	min=NULL;
	poly=NULL;
	next=NULL;
}

void Liste::insert(Liste * maillon_)
{
	maillon_->setNext(this->getNext());
	this->setNext(maillon_);
}

int Liste::compte(){
  Liste *l;
  int tmp = 0;
  l=this;
  while(l != NULL){
          tmp= tmp+1;
          l=l->getNext();
  }
  return(tmp);
}

Liste * Liste::removeDoublon()
{
  Liste *next = this->getNext();
  if(next != NULL)
  {
          if(next->getPolynome() == this->getPolynome())
          {
            //std::cerr<<"erase"<<std::endl;
                  this->setMin(next->getMin());
                  this->setNext(next->getNext());
                  next->setToNull();
                  delete next;
                  return(this);
          } else
          {
                  return(next);
          }
  } else
  {
          return(NULL);
  }
}

void Liste::checkForDoublon()
{
  Liste *l = this;
  while(l != NULL)
  {
          l=l->removeDoublon();
  }
}

void Liste::computeRoots(double a0_)
{
        Liste *l;
        l=this;
        while(l != NULL)
        {
                l->getPolynome()->roots(a0_);
                l=l->getNext();
        }
}

void Liste::add(double a2_, double a1_, double a0_)
{
        Liste *l;
        l=this;
        while(l != NULL)
        {
                l->getPolynome()->add(a2_, a1_, a0_);
                l=l->getNext();
        }
}

void Liste::computeMinOrMax(double * min, int * which)
{
        Liste *l;
        double tmp = NUMLIB_POSINF;
        *min = NUMLIB_POSINF;
        * which=-1;
        l=this;
        while(l != NULL)
        {
                l->getPolynome()->minOrMax(min, &tmp, which);
                l=l->getNext();
        }
}
#endif //INLINED
#endif
