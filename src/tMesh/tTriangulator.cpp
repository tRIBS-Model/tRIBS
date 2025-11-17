/*******************************************************************************
 * TIN-based Real-time Integrated Basin Simulator (tRIBS)
 * Distributed Hydrologic Model
 *
 * Copyright (c) 2025. tRIBS Developers
 *
 * See LICENSE file in the project root for full license information.
 ******************************************************************************/

/***************************************************************************
**
** tTriangulator.cpp  A triangulation routine based on Tipper's convex hull 
**                  algorithm Computers and geoscience vol17 no 5 pp 597-632,
**                  1991. Scaling is nlogn for random datasets. Mike Bithell 
**                  31/08/01. Arnaud Desitter. Q2 2002. Debugging and 
**                  extension of the algorithm. Binding to CHILD data 
**                  structures.
**
\***************************************************************************/

#include "src/tMesh/tTriangulator.h"



// vector product (or cross product) of p0p1,p0p2
static
double vecprod(int p0,int p1,int p2,const point *p){
	return
    (p[p1].x-p[p0].x)*(p[p2].y-p[p0].y)
    -(p[p1].y-p[p0].y)*(p[p2].x-p[p0].x);
}

const point &point::operator=( const point &p ) {
	if ( &p != this ) {
		x=p.x; y=p.y; id=p.id;
	}
	return *this;
}

// we sort according x *and* y (if x are equal) so that
// vertically aligned point are ordered.
int point::operator < (const point& p) const {
	if (x == p.x)
		return y<p.y;
	return x<p.x;
}

#if defined(DEBUG_PRINT)
void point::print () const {cout << x << ' '<< y <<endl;}
#endif
void point::write(ofstream& f) const {f<<x<<' '<<y<<endl;}


#if defined(DEBUG_PRINT)
void edge::print(const point p[]) const {p[from].print();p[to].print();}
#endif
void edge::write(ofstream& f,const point p[]) const {p[from].write(f);p[to].write(f);}
bool edge::visible(const point p[],int i) const {
	//test whether an edge on the hull is visible from a point
	//rely on the fact that a) hull is anticlockwise oriented
	//b)data is positive x ordered
	
	// that is if angle(from-to,from-i) < 0 
	const double v = vecprod(from,to,i, p);
	if (v>=0)
		return false;
	return true;
}

static
bool needswap(int i1, int i2, int i3, int i4, const point p[]){
	// i3 - i4 currently joined by diagonal.
	// i1 - i4 - i2 - i3 anti-clockwise
	// check if i1 lies inside the circumcircle for the triangle
	// i4-i2-i3 using the algorithm of Cline and Renka for roundoff
	// error.
	
	// t1=angle(i2,i3,i4), t2=angle(i1i4,i1i3), 
	// a swap is needed if t1+t2 > pi
	// Since t1+t2 < 2pi, this is equivalent to sin(t1+t2) < 0
	// This leads to: cos(t1)*cos(t2)+sin(t1)*sin(t2) < 0
	//
	// see discussion in:
	// Sloan, S.W. "A fast algorihtm for constructing Delaunay 
	// triangulations in the plane", Adv. Eng, Software,
	// 1987, 9(1)
	// Cline, A.K., Renka, R.L., "A Storage efficient method for 
	// construction of a Thiessen triangulation", Rocky Mountain 
	// Journal of Mathemetics, 1984, 14(119)
	
	const point 
    p1(p[i3]-p[i2]), p2(p[i4]-p[i2]),
    p3(p[i4]-p[i1]), p4(p[i3]-p[i1]); 
	const double cost1=p1.dot(p2), cost2=p3.dot(p4);
	
	if (cost1>=0 && cost2>=0)
		return false;
	if (cost1<0 && cost2<0)
		return true;
	const double sint1=vecprod(i2,i3,i4,p), sint2=vecprod(i1,i4,i3,p),
		sint1t2 = sint1*cost2+sint2*cost1;
	if (sint1t2 < 0)
		return true;
	return false;
}

static
int tt_swap(int tint, edge e[], const point p[]){
	
	int& to = e[tint].to;
	int& from = e[tint].from;
	int& ret = e[tint].ret;
	int& ref = e[tint].ref;
	int& let = e[tint].let;
	int& lef = e[tint].lef;
	
	//edge swapping routine - each edge has four neighbour edges
	//left and attached to from node (lef) right attached to to node (ret) etc.
	//these edges may be oriented so that their from node is that of 
	// the current edge, or not
	//this routine takes the lions share of the CPU - hull construction by comparison
	//takes much less (by about a factor of ten!)
	if (ref==-1 || lef==-1 || let==-1 || ret==-1)return 0;
	//test orientation of left and right edges - store the indices of
	//points that are not part of the current edge
	const int leftp  = (e[lef].from==from) ? e[lef].to : e[lef].from;
	const int rightp = (e[ref].from==from) ? e[ref].to : e[ref].from;
	
	if (needswap(leftp, rightp, to, from, p)){
		//now swap the left and right edges of neighbouring edges
		//taking into account orientation
		if (e[ref].from == from){
			e[ref].lef=lef;
			e[ref].let=tint;
		}else{
			e[ref].ref=tint;
			e[ref].ret=lef;
		}
		if (e[lef].from==from){
			e[lef].ref=ref;
			e[lef].ret=tint;
		}else{
			e[lef].lef=tint;
			e[lef].let=ref;
		}
		if (e[ret].to==to){
			e[ret].lef=tint;
			e[ret].let=let;
		}else{
			e[ret].ref=let;
			e[ret].ret=tint;
		}
		if (e[let].to==to){
			e[let].ref=tint;
			e[let].ret=ret;
		}else{
			e[let].lef=ret;
			e[let].let=tint;
		}
		//change the end-points for the current edge
		to=rightp;
		from=leftp;
		//re-jig the edges
		const int rf=ref;
		ref=lef;
		int rt=ret;
		ret=rf;
		const int lt=let;
		let=rt;
		lef=lt;
		//examine the neighbouring edges for delauniness recursively - this is
		//a lot more efficient than trying to swap all edges right at the end.
		tt_swap(lef,e,p);
		tt_swap(let,e,p);
		tt_swap(ref,e,p);
		tt_swap(ret,e,p);
		return 1;
	}
	return 0;
}

class cyclist{
	//a fixed size linked cyclical list using arrays
	//the code here is not yet robust necessarily to list shrinking to zero
	//number of elements
	const cyclist &operator=( const cyclist & );  // assignment operator
	cyclist( const cyclist & );
	void delNext(int list_pos);
	void add(int ej);
	void invariant() const;
	void checkRange(int) const;
public:
		cyclist(int s);
	~cyclist();
	int getEdge(int list_pos) const {
		checkRange(list_pos);
		return ejs[list_pos].data;
	}
	int delNextPos(int list_pos);
	int delNextNeg(int list_pos);
	int getNextPos(int list_pos) const {
		checkRange(list_pos);
		return ejs[list_pos].next;
	}
	int getNextNeg(int list_pos) const{
		checkRange(list_pos);
		return ejs[list_pos].prev;
	}
	int addBefore(int a /*location*/, int ej /*data*/);
	int addAfter(int a /*location*/, int ej /*data*/);
#if defined(DEBUG_PRINT)
	void print() const;
#endif
private:
		struct item{
			int next,prev,data;
		};
	item* ejs;  // array of size "size"
	const int size; // maximum size
					// "hole" points to the next unfilled location,
					// "prev" to the location that was last filled
					// we keep track of empty bits of the list using hole to point to an empty slot,
					// and the value ejs[hole] to point to the next empty slot
					// "num" is the current number of elements 
	int prev, hole, num;
};

cyclist::cyclist(int s): ejs(0),size(s),prev(0),hole(0),num(0) {
	ejs=new item[size];
	//fill ej array with a set of pointers to the next unfilled location
	for (int i=0;i<size-1;i++) ejs[i].data=i+1;
	ejs[size-1].data = -1;
}
cyclist::~cyclist(){
	invariant();
	delete [] ejs;
}
void cyclist::checkRange(int p) const {
	assert(p<size);
	assert(0<=p);
}
void cyclist::invariant() const {
	assert(num<=size);
	assert(0<=num);
	checkRange(prev);
	assert(hole<size);
	if (num!=size)
		assert(0<=hole);
}
void cyclist::delNext(int list_pos){
	checkRange(list_pos);
	assert(num!=0);
	ejs[list_pos].data=hole;
	hole=list_pos;
	ejs[ejs[list_pos].prev].next=ejs[list_pos].next;
	ejs[ejs[list_pos].next].prev=ejs[list_pos].prev;
	num--;
}
int cyclist::delNextPos(int list_pos){
	delNext(list_pos);
	return ejs[list_pos].next;
}
int cyclist::delNextNeg(int list_pos){
	delNext(list_pos);
	return ejs[list_pos].prev;
}
void cyclist::add(int ej){
	//build hull from scratch in numerical order - we assume you got the orientation
	//right!! (anti-clockwise)
	checkRange(hole<size);
	const int n=ejs[hole].data;
	ejs[hole].data=ej;
	//prev stores the location of the place in the array that was most recently filled
	ejs[prev].next=hole;
	//rev is the set of backward pointers
	ejs[hole].prev=prev;
	//the list is cyclic
	ejs[hole].next=0;
	ejs[0].prev=hole;
	prev=hole;
	hole=n;
	num++;
}
int cyclist::addBefore(int a, int ej){
	checkRange(a);
	//first check for the empty list
	if (num ==0){add(ej);return prev;}
	//otherwise add on before the specified position, using the empty storage slot
	const int n=ejs[hole].data;
	ejs[hole].prev=ejs[a].prev;
	ejs[hole].next=a;
	ejs[ejs[a].prev].next=hole;
	ejs[a].prev=hole;
	ejs[hole].data=ej;
	prev=hole;
	hole=n;
	num++;
	//return the value that hole had at the start of the method 
	return prev;
}
int cyclist::addAfter(int a,int ej){
	checkRange(a);
	//first check for the empty list
	if (num ==0){add(ej);return prev;}
	//otherwise add on after the specified position, using the empty storage slot
	const int n=ejs[hole].data;
	ejs[hole].next=ejs[a].next;
	ejs[ejs[a].next].prev=hole;
	ejs[a].next=hole;
	ejs[hole].prev=a;
	ejs[hole].data=ej;
	prev=hole;
	hole=n;
	num++;
	//return the value that hole had at the start of the method 
	return prev;
}
#if defined(DEBUG_PRINT)
void cyclist::print() const {
	int j=ejs[0].next;
	for (int i=0;i<num;i++){
		cout<<ejs[j].data<<endl;
		j=ejs[j].next;
	}
}
#endif

// if the first points are aligned, then we build "by hand" the triangulation.
// Draw a picture to find out. Not that difficult.
void start_aligned_point(int &lower_hull_pos, int &upper_hull_pos, int &next_edge,
						 int &next_point,
						 int npoints, const point *p,  edge *edges,
						 cyclist &hull){
	// the 3 first points are aligned.
	// find the first non aligned point
	int orient = 0;
	int j;
	for(j=3;j<npoints;++j){
		const double v = vecprod(j-2,j-1,j,p);
		if (v!=0.){
			orient = -1;
			if (v>0.) 
				orient = 1;
			break;
		}
	}
	if (orient == 0){
		cout << "All nodes aligned. Ill-conditioned problem. Bailing out." 
		<< endl;
		tt_error_handler();
	}
	// if orient > 0, build edge from 0 to j-1
	// if orient < 0, build edge from j-1 to 0
	//if (0) // DEBUG
	//  cout << "first non aligned j=" << j << " o=" << orient << endl;
	
	const int nap = j; // first non aligned point
	if (orient > 0) {
		// aligned edges
		int inode = 0;
		for(int iedge = 0; iedge<nap-1; ++iedge ){
			edges[iedge].from=inode;
			edges[iedge].to=inode+1;
			inode++;
			//if (0) // DEBUG
			//cout << "edge=" << iedge << " from=" << edges[iedge].from << " to=" 
			//<< edges[iedge].to << endl; 
		}
		
	} else {
		int inode = nap-2;
		for(int iedge = 0; iedge<nap-1; ++iedge ){
			edges[iedge].from=inode+1;
			edges[iedge].to=inode;
			inode--;
			//if (0) // DEBUG
			//cout << "edge=" << iedge << " from=" << edges[iedge].from << " to=" 
			//<< edges[iedge].to << endl; 
		}
	}
	// close the domain
	edges[nap-1].from=edges[nap-2].to;
	edges[nap-1].to=nap;
	edges[nap].from=nap;
	edges[nap].to=edges[0].from;
	// fill cyclic list
	{
		// because the triangulation done, lower_hull_pos and upper_hull_pos
		// should be properly defined
		int start=hull.addAfter(0,0);
		for (int iedge=1;iedge<nap-1; ++iedge){
			start=hull.addAfter(start,iedge);
		}
		start=hull.addAfter(start,nap-1);
		lower_hull_pos=start;
		start=hull.addAfter(start,nap);
		upper_hull_pos=start;
	}
	// add interior edges (nap-2 of them)
	{
		int iedge2=0;
		for(int iedge=nap+1;iedge<2*nap-1;++iedge){
			edges[iedge].from=edges[iedge2].to;
			edges[iedge].to=nap;
			iedge2++;
		}
	}
	// compute neighbour edges
	for(int iedge=0;iedge<nap-2;++iedge){
		edges[iedge].lef=nap+iedge;
		edges[iedge].let=nap+1+iedge;
	}
	edges[nap-2].lef=2*nap-2;
	edges[nap-2].let=nap-1;
	// closure of domain
	edges[nap-1].lef=nap-2;
	edges[nap-1].let=2*nap-2;
	edges[nap].lef=nap+1;
	edges[nap].let=0;
	// interior edges
	edges[nap+1].lef=0;
	edges[nap+1].let=nap;
	edges[nap+1].ref=1;
	if (nap!=3){ // more than one interior edge
		edges[nap+1].ret=nap+2;
		int iedge2 = 1;
		for(int iedge=nap+2;iedge<2*nap-2;++iedge){
			edges[iedge].lef=iedge2;
			edges[iedge].let=iedge-1;
			edges[iedge].ref=iedge2+1;
			edges[iedge].ret=iedge+1;
			iedge2++;
		}
		edges[2*nap-2].lef=iedge2;
		edges[2*nap-2].let=2*nap-3;
		edges[2*nap-2].ref=iedge2+1;
	}
	edges[2*nap-2].ret=nap-1;
	//
	next_edge = 2*nap-1; // number of existing edges: 2*nap-1
	next_point = nap+1;
	//if (0) // DEBUG
	//cout << "next point=" << next_point << endl;
	//
	//if (0) { //DEBUG
	// for (int iedge=0;iedge<next_edge; iedge++){
	//   cout << "iedge=" << iedge
	//	   << " from=" << edges[iedge].from
	//	   << " to=" << edges[iedge].to
	//	   << " lef=" << edges[iedge].lef 
	//	   << " let=" << edges[iedge].let 
	//	   << " ref=" << edges[iedge].ref 
	//	   << " ret=" << edges[iedge].ret 
	//	   << endl;
	//}
	//}
}

void triangulate(int npoints,const point p[], int *pnedges, edge** edges_ret){
	
	assert(npoints >= 3);
	
	//convex hull is a cyclical list - it will consist of anticlockwise
	//ordered edges - since each new point adds at most 1 extra edge (nett)
	//to the hull, there are at most npoints edges on the hull
	cyclist hull(npoints);
	
	//and the edges - there are at most three edges per point
	const long nn=3*npoints;
	edge* edges;
	edges=new edge[nn];
	
	int lower_hull_pos, upper_hull_pos, next_edge, next_point;
	
	//make first three edges  - these will form the initial convex hull
	//make sure orientation is anticlockwise
	// Arnaud: orientation is counter clockwise <=> angle(p0p1,p0p2) >= 0
	//  <=> sin(p0p1,p0p2) >=0 <=> vect_prod(p0p1,p0p2) >= 0
	{
		const double v = vecprod(0,1,2,p);
		if (v!=0.) {
			if (v>0.) {
				edges[0].from=0;
				edges[0].to=1;
				edges[1].from=1;
				edges[1].to=2;
				edges[2].from=2;
				edges[2].to=0;
			} else {
				edges[0].from=0;
				edges[0].to=2;
				edges[1].from=2;
				edges[1].to=1;
				edges[2].from=1;
				edges[2].to=0;
			}
			//make left edges
			edges[0].lef=2;
			edges[0].let=1;
			edges[1].lef=0;
			edges[1].let=2;
			edges[2].lef=1;
			edges[2].let=0;
			// add the edges to the hull in order
			// get upper and lower edges as indices into the hull
			// 'cos hull is cyclic, we don't need to bother if these are visible
			// but they do need to be oriented so that upper edge is further round the 
			// hull in a positive direction than lower_hull_pos
			int start=hull.addAfter(0,0);
			lower_hull_pos=hull.addAfter(start,1); 
			upper_hull_pos=hull.addAfter(lower_hull_pos,2);
			// loop through the remaining points adding edges to the hull
			next_edge=3;
			next_point=3;
		} else {
			start_aligned_point(lower_hull_pos, upper_hull_pos, next_edge,
								next_point,
								npoints, p, edges, hull);
		}
	}
	
	int count,saved_edge;
	
	for (int i=next_point;i<npoints;i++){
		saved_edge=-1;
		if (p[i-1].x == p[i].x && p[i-1].y == p[i].y){
			continue; // skip the duplicated point
		}
		
		//go round the hull looking for visible edges - we need to go round
		//in two directions from the current upper and lower edges
		//first set up the new edge that joins to the point coincident
		//between upper and lower edge
		if(edges[hull.getEdge(upper_hull_pos)].visible(p,i)){
			const int hup = hull.getEdge(upper_hull_pos);
			//make new edge - from and to nodes preserve hull orientation
			edges[next_edge].from=edges[hup].from;
			edges[next_edge].to=i;
			//save hull position for possible later use
			saved_edge=next_edge;
			//connectivity for swapping - we know the ID of the
			//next edge to be created since upper hull is visible
			edges[next_edge].lef=hup;
			edges[next_edge].let=next_edge+1;
			next_edge++;
		}
		else{
			//we can't see the upper edge - can we see the lower one?
			if(!edges[hull.getEdge(lower_hull_pos)].visible(p,i)){
				//can't see the upper or lower edge - chose a bad initial state! go round
				//the hull a bit.Upper hull will still be invisible (its what used to be lower hull)
				if (i==3){ // only when starting from a single initial triangle
					lower_hull_pos=hull.getNextPos(upper_hull_pos);
					upper_hull_pos=hull.getNextPos(lower_hull_pos);
				}else{
					//or else its an error!
					cout<<"Triangulate: Can't see the hull from the new point!? number is "<<i<<endl;
					tt_error_handler();
				}
			}
			const int hlow = hull.getEdge(lower_hull_pos);
			edges[next_edge].from=i;
			edges[next_edge].to=edges[hlow].to;
			
			//connectivity for swapping - we know the ID of the
			//next edge but one to be created since upper hull is *not* visible
			edges[next_edge].let=hlow;
			edges[next_edge].lef=next_edge+1;
			next_edge++;
		}
		//now we need to add the upper hull edges - drop through to set edge made above
		//to the new upper edge if upper edge not visible
		count=0;
		while (edges[hull.getEdge(upper_hull_pos)].visible(p,i)){
			const int hup=hull.getEdge(upper_hull_pos);
			edges[next_edge].from=i;
			edges[next_edge].to=edges[hup].to;
			//if we got here, the upper hull is visible, so we know which way the edges point
			//set the left and right neighbour edges
			if (count!=0){
				//if we're on the second time round the loop, the edge made in
				// the previous pass needs to be connected on its right side
				edges[next_edge-1].ref=next_edge;
				edges[next_edge-1].ret=hup;
			}
			count++;
			edges[hup].ref=next_edge-1;
			edges[hup].ret=next_edge;
			edges[next_edge].let=hup;
			edges[next_edge].lef=next_edge-1;
			//check the hull edge's delauniness
			tt_swap(hup,edges,p);
			next_edge++;
			//delete upper edge from the hull
			//and go round the hull in the positive direction
			upper_hull_pos=hull.delNextPos(upper_hull_pos);
		}
		//we drop through to here if upper edge is not visible
		//add edge to the hull
		//set the new upper edge to the most recently created upper edge
		upper_hull_pos=hull.addBefore(upper_hull_pos,next_edge-1);
		//now we need to add the lower hull edges
		while (edges[hull.getEdge(lower_hull_pos)].visible(p,i)){
			const int hlow=hull.getEdge(lower_hull_pos);
			//upper hull was not visible
			if (saved_edge==-1){
				//after this pass through the loop, always do the else clause
				saved_edge=next_edge-1;
			}else{
				//connect right side of previously created edge
				edges[saved_edge].ref=hlow;
				edges[saved_edge].ret=next_edge;
			}
			edges[next_edge].to=i;
			edges[next_edge].from=edges[hlow].from;
			edges[hlow].ret=saved_edge;
			edges[hlow].ref=next_edge;
			edges[next_edge].lef=hlow;
			edges[next_edge].let=saved_edge;
			//swap it if it needs it
			tt_swap(hlow,edges,p);
			//keep the edge for below - this is necessary in case the upper hull wasn't visible
			saved_edge=next_edge;
			next_edge++;
			//delete lower edge from the hull
			//and go round the hull in the negative direction
			lower_hull_pos=hull.delNextNeg(lower_hull_pos);
		}
		//add edge to the hull
		//set the new lower edge to the most recently created lower edge
		//if no lower hull pos was visible, use the saved edge from the upper pos earlier on
		lower_hull_pos=hull.addAfter(lower_hull_pos,saved_edge);
	}
	
	// results
	{ 
		int i=0;
		while(edges[i].from != -1 ) i++;
		*pnedges = i;
	}
	*edges_ret = edges;
	
}

#include "heapsort.h"

static
void tt_verify_sort(int npoints, const point *p){
	for(int ipoint=1; ipoint < npoints; ++ipoint){
		assert(p[ipoint-1].x <= p[ipoint].x);
		if (p[ipoint-1].x == p[ipoint].x)
			assert(p[ipoint-1].y <= p[ipoint].y);
	}
}

void tt_sort_triangulate(int npoints, point *p,
						 int *pnedges, edge** edges_ret){
	
	//sort the points - note that the point class defines the
	// < operator so that the sort is on the x co-ordinate
	//array p will be replaced with the array sorted in x
	heapsort(npoints,p);
	if (0)
		tt_verify_sort(npoints, p);
	
#if defined(TIMING)
	{
		time_t t1 = time(NULL);
		clock_t tick1 = clock();
#endif
		
		//triangulate the set of points
		triangulate(npoints,p,pnedges, edges_ret);
		
		
#if defined(TIMING)
		time_t t2 = time(NULL);
		clock_t tick2 = clock();
		cout << "elapsed time (time)= " << difftime(t2,t1) << " s"
			<< " (clock)= " << static_cast<double>(tick2-tick1)/CLOCKS_PER_SEC << " s"
			<< endl;
	}
#endif
	/*if (0) { // DEBUG
		for (int iedge=0;iedge<*pnedges; iedge++){
			cout << "iedge=" << iedge
			<< " from=" << (*edges_ret)[iedge].from
			<< " to=" << (*edges_ret)[iedge].to
			<< " lef=" << (*edges_ret)[iedge].lef 
			<< " let=" << (*edges_ret)[iedge].let 
			<< " ref=" << (*edges_ret)[iedge].ref 
			<< " ret=" << (*edges_ret)[iedge].ret 
			<< endl;
		}
	}*/
}

// Auxilary class used when building the element to node connectivity
// table
class edge_auxi_t {
public:
	edge_auxi_t()
    :
    left_visited_(false), right_visited_(false),
    ie_left(-2), ie_right(-2)
	{}
	bool left_visited() const { return left_visited_; }
	bool right_visited() const { return right_visited_; }
	int ielem_left() const;
	int ielem_right() const;
	void mark_left(int ielem);
	void mark_right(int ielem);
private:
		bool left_visited_, right_visited_;
	// element on the right and left side
	int ie_left, ie_right;
};

void edge_auxi_t::mark_left(int ielem) { 
	left_visited_ = true;
	ie_left = ielem;
}
void edge_auxi_t::mark_right(int ielem) { 
	right_visited_ = true;
	ie_right = ielem;
} 
int edge_auxi_t::ielem_left() const { assert(left_visited_); return ie_left; }
int edge_auxi_t::ielem_right() const { assert(right_visited_); return ie_right; }


const oriented_edge& oriented_edge::operator=( const oriented_edge &_e ){
	if (this != &_e) {
		_edge = _e.e(); _orientation = _e.o();
	}
	return *this;
}

oriented_edge oriented_edge::next_ccw_around_from(const edge* edges) const {
	int ires;
	bool bres = true;
	if (o()){
		ires = edges[e()].lef;
	} else {
		ires = edges[e()].ret;
	}
	if (ires != -1){
		bres = (edges[ires].ref == e() ? true:false);
		if (!bres)
			assert( edges[ires].let == e());
	}
	return oriented_edge(ires,bres);
}

oriented_edge oriented_edge::next_cw_around_from(const edge* edges) const {
	int ires;
	bool bres = true;
	if (o()){
		ires = edges[e()].ref;
	} else {
		ires = edges[e()].let;
	}
	if (ires != -1){
		bres = (edges[ires].lef == e() ? true:false);
		if (!bres)
			assert( edges[ires].ret == e());
	}
	return oriented_edge(ires,bres);
}

// give counter clockwise edge
oriented_edge oriented_edge::ccw_edge_around_from(const edge* edges) const {
	oriented_edge ccw = next_ccw_around_from(edges);
	if (ccw.nonvalid()) {
		// iterated on clockwise edges
		oriented_edge e1(*this);
		for(;;){
			const oriented_edge enext = e1.next_cw_around_from(edges);
			if (enext.nonvalid()) 
				break;
			e1 = enext;
		}
		ccw = e1;
	}
	return ccw;
}

// mark as visited the side of iedge_markable that points to
// iedge_orig
static
void mark_as_visited(int iedge_markable, int iedge_orig,
					 const edge* edges, edge_auxi_t* edges_visit,
					 int ielem){
	if (edges[iedge_markable].lef == iedge_orig ||
		edges[iedge_markable].let == iedge_orig){
		edges_visit[iedge_markable].mark_left(ielem);
	} else {
		assert(edges[iedge_markable].ref == iedge_orig ||
			   edges[iedge_markable].ret == iedge_orig);
		edges_visit[iedge_markable].mark_right(ielem);
	}
}

// check geometry properties of the edges in the elems table
static
void sanity_check_elems(const edge* edges, const elem *elems, 
						const edge_auxi_t *edges_visit, int ielem){
	if (elems[ielem].eo2){
		assert(edges[elems[ielem].e2].from == elems[ielem].p2);
		assert(edges[elems[ielem].e2].to == elems[ielem].p1);
		assert(edges_visit[elems[ielem].e2].ielem_right() == ielem);
	} else {
		assert(edges[elems[ielem].e2].from == elems[ielem].p1);
		assert(edges[elems[ielem].e2].to == elems[ielem].p2);
		assert(edges_visit[elems[ielem].e2].ielem_left() == ielem);
	}
	if (elems[ielem].eo1){
		assert(edges[elems[ielem].e1].from == elems[ielem].p1);
		assert(edges[elems[ielem].e1].to == elems[ielem].p3);
		assert(edges_visit[elems[ielem].e1].ielem_right() == ielem);
	} else {
		assert(edges[elems[ielem].e1].from == elems[ielem].p3);
		assert(edges[elems[ielem].e1].to == elems[ielem].p1);
		assert(edges_visit[elems[ielem].e1].ielem_left() == ielem);
	}
	if (elems[ielem].eo3){
		assert(edges[elems[ielem].e3].from == elems[ielem].p3);
		assert(edges[elems[ielem].e3].to == elems[ielem].p2);
		assert(edges_visit[elems[ielem].e3].ielem_right() == ielem);
	} else {
		assert(edges[elems[ielem].e3].from == elems[ielem].p2);
		assert(edges[elems[ielem].e3].to == elems[ielem].p3);
		assert(edges_visit[elems[ielem].e3].ielem_left() == ielem);
	}
}

// build connectivity table element to node and edges
// build pelem, pelems_ret
void tt_build_elem_table(int npoints, const point *p,
						 int nedges, const edge* edges,
						 int *pnelem, elem** pelems_ret){
	// Euler invariant for a triangulation:
	const int nelem = 1 + nedges - npoints;
	elem *elems = new elem[nelem]; 
	edge_auxi_t *edges_visit = new edge_auxi_t[nedges];
	{ 
		// build non oriented edges for each element
		// scan all edges. Assign element to right and left
		// and mark the visited sides for each edges
		
		// build edges per element
		int ielem = 0;
		for(int iedge=0;iedge<nedges;iedge++) {
			// left
			if (! edges_visit[iedge].left_visited()) {
				int ielem_current = -1;
				if (edges[iedge].lef == -1) {
					assert(edges[iedge].let == -1);
				}
                else {
					// don't bother with orientation at the moment
					ielem_current = ielem;
					elems[ielem].e1 = iedge;
					elems[ielem].e2 = edges[iedge].lef;
					mark_as_visited(edges[iedge].lef,iedge,edges, edges_visit, ielem_current);
					elems[ielem].e3 = edges[iedge].let;
					mark_as_visited(edges[iedge].let,iedge,edges, edges_visit, ielem_current);
					ielem++;
				}
				edges_visit[iedge].mark_left(ielem_current);
			}
			// right
			if (! edges_visit[iedge].right_visited()) {
				int ielem_current = -1;
				if (edges[iedge].ref == -1) {
					assert(edges[iedge].ref == -1);
				} else {
					// don't bother with orientation at the moment
					ielem_current = ielem;
					elems[ielem].e1 = iedge;
					elems[ielem].e2 = edges[iedge].ref;
					mark_as_visited(edges[iedge].ref,iedge,edges, edges_visit, ielem_current);
					elems[ielem].e3 = edges[iedge].ret;
					mark_as_visited(edges[iedge].ret,iedge,edges, edges_visit, ielem_current);
					ielem++;
				}
				edges_visit[iedge].mark_right(ielem_current);
			}
		}
		assert(ielem == nelem);
	}
	for(int ielem=0;ielem<nelem;ielem++){
    {
		// build vertices per element
#define SWAP_E(EI,EJ) \
		do { int etemp = elems[ielem].EI; \
			elems[ielem].EI = elems[ielem].EJ; \
				elems[ielem].EJ = etemp; } while(0)
					
					elems[ielem].p1 = edges[elems[ielem].e2].to;
		elems[ielem].p2 = edges[elems[ielem].e2].from;
		// using e3
		int itemp = edges[elems[ielem].e3].from;
		if (itemp == elems[ielem].p1 || itemp == elems[ielem].p2){
			elems[ielem].p3 = edges[elems[ielem].e3].to;
			if (itemp == elems[ielem].p1) { // e3 joins p1-p3
				SWAP_E(e1,e3);
			}
		} else {
			elems[ielem].p3 = itemp;
			if (edges[elems[ielem].e3].to == elems[ielem].p1){ // e3 joins p1-p3
				SWAP_E(e1,e3);
			}
		}
		// orientation: counter clockwise
		// angle(p2p1.p2p3) must be negative (counter clockwise)
		// <=> sin(p2p1.p2p3) <= 0 <=> vect_prod(p2p1,p2p3) <=0
		const double v=
			vecprod(elems[ielem].p2,elems[ielem].p1,elems[ielem].p3,p);
		// points should not be aligned.
		if (v==0.){
			cout << "These points are aligned: "
			<< elems[ielem].p1 << "(" << p[elems[ielem].p1].x
			<< "," << p[elems[ielem].p1].y << "), " 
			<< elems[ielem].p2 << "(" << p[elems[ielem].p2].x
			<< "," << p[elems[ielem].p2].y << "), " 
			<< elems[ielem].p3 << "(" << p[elems[ielem].p3].x
			<< "," << p[elems[ielem].p3].y << ")" 
			<< endl;
		}
		assert(v != 0.);
		if (v>0) {
			SWAP_E(p1,p3);
			SWAP_E(e2,e3);
		}
		// edges orientations
		elems[ielem].eo2 = ((edges[elems[ielem].e2].from == elems[ielem].p2) ? true:false);
		elems[ielem].eo3 = ((edges[elems[ielem].e3].from == elems[ielem].p3) ? true:false);
		elems[ielem].eo1 = ((edges[elems[ielem].e1].from == elems[ielem].p1) ? true:false);
#undef SWAP_E
    }
		if (1) //DEBUG
			sanity_check_elems(edges, elems, edges_visit, ielem);
		{
			// build t1, t2, t3
			if (elems[ielem].eo2){
				elems[ielem].t3 = edges_visit[elems[ielem].e2].ielem_left();
			} else {
				elems[ielem].t3 = edges_visit[elems[ielem].e2].ielem_right();
			}
			if (elems[ielem].eo1){
				elems[ielem].t2 = edges_visit[elems[ielem].e1].ielem_left();
			} else {
				elems[ielem].t2 = edges_visit[elems[ielem].e1].ielem_right();
			}
			if (elems[ielem].eo3){
				elems[ielem].t1 = edges_visit[elems[ielem].e3].ielem_left();
			} else {
				elems[ielem].t1 = edges_visit[elems[ielem].e3].ielem_right();
			}
			
		}
	}
	/*if (0) { //DEBUG
		for(int ielem=0;ielem<nelem;ielem++){
			cout << "elem=" << ielem
			<< " p1=" << elems[ielem].p1
			<< " (" << p[elems[ielem].p1].x <<"," << p[elems[ielem].p1].y << ")"
			<< " p2=" << elems[ielem].p2
			<< " (" << p[elems[ielem].p2].x <<"," << p[elems[ielem].p2].y << ")"
			<< " p3=" << elems[ielem].p3
			<< " (" << p[elems[ielem].p3].x <<"," << p[elems[ielem].p3].y << ")"
			<< endl;
		}
    }*/
	delete [] edges_visit; edges_visit = NULL;
	*pnelem = nelem;
	*pelems_ret = elems;
}

void tt_build_spoke(int npoints, int nedges, const edge* edges,
					oriented_edge** poedge){
	*poedge = new oriented_edge[npoints];
	for(int iedge=0;iedge!=nedges;++iedge){
		const int from=edges[iedge].from;
		if ((*poedge)[from].nonvalid()) {
			(*poedge)[from].set(iedge, true);
		} else {
			const int to=edges[iedge].to;
			if ((*poedge)[to].nonvalid()) {
				(*poedge)[to].set(iedge, false);
			}
		}
	}
}

void tt_sort_triangulate(int npoints, point *p,
						 int *pnedges, edge** edges_ret,
						 int *pnelem, elem** pelems_ret){
	tt_sort_triangulate(npoints,p,pnedges, edges_ret);
	tt_build_elem_table(npoints, p, *pnedges, *edges_ret, pnelem, pelems_ret);
}


void tt_error_handler(void){
	cout<<"\nFatal Error in tTriangulator";
}


//=========================================================================
//
//
//                       End of tTriangulator.cpp
//
//
//=========================================================================
