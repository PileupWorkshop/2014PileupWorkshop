#include "VoronoiParticleArea.hh"

#ifndef NOCGAL
// Delaunay triangulation (as in FJ)
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Voronoi graph
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#endif //NOCGAL

using namespace std;
using namespace fastjet;


#ifndef NOCGAL
// cgal info (copied from FJ)
//------------------------------------------------------------------------
/// A class to provide an "int" with an initial value.
const int NEW_VERTEX=-1;
class InitialisedInt {
 private:
  int _val;
 public:
  inline InitialisedInt () {_val=NEW_VERTEX;};
  inline InitialisedInt& operator= (int value) {_val = value; return *this;};
  inline int val() const {return _val;};
};


struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Triangulation_vertex_base_with_info_2<InitialisedInt,K> Vbb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb> Vb;
typedef CGAL::Triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>  Dt;
typedef CGAL::Triangulation_hierarchy_2<Dt> Triangulation;
#endif //NOCGAL

void VoronoiParticleAreas::compute_areas(const vector<PseudoJet> &event){
#ifndef NOCGAL
  Triangulation TR;
  vector<Triangulation::Vertex_handle> vertices(3*event.size());
                       
  // create the original voronoi graph
  Triangulation::Vertex_handle v_copy1, v_copy2;
  for (unsigned int i=0;i<event.size();i++){
    double rap = event[i].rap();
    double phi = event[i].phi();
    Triangulation::Vertex_handle v = TR.insert(Triangulation::Point(rap, phi));
    if (v->info().val() != NEW_VERTEX) {
      ostringstream err;
      cerr << "CGAL Warning: point "<< i <<" coincides with point "
           << v->info().val() << " discarding the copy" << endl;
    } else {
      v->info() = i;
      
      // insert copies to simulate periodicity in y
      v_copy1 = TR.insert(Triangulation::Point(rap, phi-twopi));
      v_copy1->info() = i;
      v_copy2 = TR.insert(Triangulation::Point(rap, phi+twopi));
      v_copy2->info() = i;
      
    }
    vertices[i] = v;
    vertices[i+event.size()] = v_copy1;
    vertices[i+2*event.size()] = v_copy2;
  }

  // loop over vertices and compute the areas
  K::Iso_rectangle_2 acceptance(-ymax, -4.0*pi, ymax, 6.0*pi);
  vector<K::Point_2> intersections;

  areas.clear();
  areas.resize(event.size());

  for (unsigned int i=0; i< event.size(); i++){
    intersections.clear();

    // circulate over the edges and build the cell intersections
    Triangulation::Edge_circulator ecirc0 = TR.incident_edges(vertices[i]);
    Triangulation::Edge_circulator ecirc = ecirc0;
    do{
      // get the Voronoi dual
      if (! (TR.is_infinite(ecirc))){
        CGAL::Object vd_edge = TR.dual(ecirc);
        CGAL::Object edge_intersection;
          
        // check what type of connection this is
        //  (i) a segment
        const K::Segment_2* seg=CGAL::object_cast<K::Segment_2>(&vd_edge);
        if (seg){
          edge_intersection = CGAL::intersection(*seg, acceptance);
        } else {
          //  (ii) a ray (half-line)
          const K::Ray_2* ray=CGAL::object_cast<K::Ray_2>(&vd_edge);
          if (ray){
            edge_intersection = CGAL::intersection(*ray, acceptance);
          } else {
            //  (iii) a line (unbounded)
            const K::Line_2* line=CGAL::object_cast<K::Line_2>(&vd_edge);
            if (line)
              edge_intersection = CGAL::intersection(*line, acceptance);
          }
        }
      
        // if we have an intersection, keep the vertices for future area computation
        const K::Segment_2 *seg_inside = CGAL::object_cast<K::Segment_2>(&edge_intersection);
        if (seg_inside){
          intersections.push_back(seg_inside->vertex(0));
          intersections.push_back(seg_inside->vertex(1));
        }
      }
      ecirc++;
    } while (ecirc != ecirc0);

    // loop over all the intersections and compute the area
    double area=0.0;
    unsigned int n_intersections = intersections.size();
    assert(n_intersections>2); // should at least be a triangle
    const K::Point_2 & centre = vertices[i]->point();
    for (unsigned int j=0; j< n_intersections-1; j++)
      area+=K::Triangle_2(centre, intersections[j], intersections[j+1]).area();
    area+=K::Triangle_2(centre, intersections[n_intersections-1], intersections[0]).area();

    areas[i] = area;
  }
#endif // NOCGAL
}


