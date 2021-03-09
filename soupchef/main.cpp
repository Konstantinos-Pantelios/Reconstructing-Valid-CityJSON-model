#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include<map>
#include <unordered_map>
#include "DCEL.hpp"


#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmultichar"
// forward declarations; these functions are given below main()
void DemoDCEL();
void printDCEL(DCEL & D);

//~~~~~~~~~~~~~~~~~~~~~ 09-03-2021 Read .obj file into memory - vertices and faces ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//import from hw1 code
void read(const char *file_in, std::vector<std::vector<double>> &v, std::vector<std::vector<unsigned int>> &f) {
    std::ifstream stream_in;
    stream_in.open(file_in);
    if (stream_in.is_open()) {
        std::string line;
        while (getline(stream_in, line)) {
            std::istringstream iss(line);
            std::string word;
            iss >> word;
            if (word == "v") {
                std::vector<double> coordinates;
                while (iss >> word)
                    coordinates.push_back(std::stof(word));
                if (coordinates.size() == 3) v.push_back({coordinates[0], coordinates[1], coordinates[2]});
            } else if (word == "f") {
                f.push_back(std::vector<unsigned int>());
                while (iss >> word) f.back().push_back((unsigned int) std::stoul(word));
            }
        }
    }
}
//~~~~~~~~~~~~~~~~~~~~~ 09-03-2021 Read .obj file into memory - vertices and faces ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/* 
  Example functions that you could implement. But you are 
  free to organise/modify the code however you want.
  After each function you should have a DCEL without invalid elements!
*/
// 1.
void importOBJ(DCEL & D, const char *file_in) {
  std::vector<std::vector<double>> vertices;
  std::vector<std::vector<unsigned int>> faces;
    read(file_in,vertices,faces);
  /*
  std::unordered_map<unsigned int, Vertex*> hashmap;
  read(file_in,vertices,faces);
  unsigned int c = 0;
  for(auto const& i : vertices){
      Vertex* v = D.createVertex(i[0],i[1],i[2],c);
      hashmap.insert(std::make_pair(c,v));
      c++;
  }
  for(auto const& j: faces){
      Face* f = D.createFace();
      HalfEdge* e0 = D.createHalfEdge();
      HalfEdge* e1 = D.createHalfEdge();
      HalfEdge* e2 = D.createHalfEdge();
      e0->origin = hashmap[j[0]-1];
      e0->destination = hashmap[j[1]-1];
      e0->incidentFace = f; */
    for (unsigned int i = 0; i < faces.size(); i++) {
        Vertex* v0 = D.createVertex(vertices[(faces[i][0])-1][0], vertices[(faces[i][0])-1][1], vertices[(faces[i][0])-1][2], faces[i][0]-1); //check indices for 0-1
        Vertex* v1 = D.createVertex(vertices[(faces[i][1])-1][0], vertices[(faces[i][1])-1][1], vertices[(faces[i][1])-1][2], faces[i][1]-1);
        Vertex* v2 = D.createVertex(vertices[(faces[i][2])-1][0], vertices[(faces[i][2])-1][1], vertices[(faces[i][2])-1][2], faces[i][2]-1);

        HalfEdge* e0 = D.createHalfEdge();
        HalfEdge* e1 = D.createHalfEdge();
        HalfEdge* e2 = D.createHalfEdge();
        HalfEdge* e3 = D.createHalfEdge();
        HalfEdge* e4 = D.createHalfEdge();
        HalfEdge* e5 = D.createHalfEdge();

        Face* fi = D.createFace();

        e0->origin = v0;
        e0->destination = v1;
        e0->twin = e3;
        e0->next = e1;
        e0->prev = e2;
        e0->incidentFace = fi;

        e3->origin = v1;
        e3->destination = v0;
        e3->twin = e0;
        e3->next = e5;
        e3->prev = e4;


        //If a half-edge is incident to 'open space' (ie not an actual face with an exterior boundary),
        //we use the infiniteFace which is predifined in the DCEL class

        e3->incidentFace = D.infiniteFace();

        e1->origin = v1;
        e1->destination = v2;
        e1->twin = e4;
        e1->next = e2;
        e1->prev = e0;
        e1->incidentFace = fi;

        e4->origin = v2;
        e4->destination = v1;
        e4->twin = e1;
        e4->next = e3;
        e4->prev = e5;
        e4->incidentFace = D.infiniteFace();

        e2->origin = v2;
        e2->destination = v0;
        e2->twin = e5;
        e2->next = e0;
        e2->prev = e1;
        e2->incidentFace = fi;

        e5->origin = v0;
        e5->destination = v2;
        e5->twin = e2;
        e5->next = e4;
        e5->prev = e3;
        e5->incidentFace = D.infiniteFace();

        fi->exteriorEdge = e0;
    }

  }


// 2.
void groupTriangles(DCEL & D) {
  // to do
}
// 3.
void orientMeshes(DCEL & D) {
  // to do
}
// 4.
void mergeCoPlanarFaces(DCEL & D) {
  // to do
}
// 5.
void exportCityJSON(DCEL & D, const char *file_out) {
  // to do
}


int main(int argc, const char * argv[]) {
  const char *file_in = "/home/konstantinos/Desktop/TUDelft-Courses/Q3/GEO1004/hw2/cube_soup.obj";
  const char *file_out = "/home/konstantinos/Desktop/TUDelft-Courses/Q3/GEO1004/hw2/cube.json";




  // Demonstrate how to use the DCEL to get you started (see function implementation below)
  // you can remove this from the final code
  //DemoDCEL();

  // create an empty DCEL
  DCEL D;
  // 1. read the triangle soup from the OBJ input file and convert it to the DCEL,

    importOBJ(D,file_in);
    D.vertices();
  //~~~~~~~~~~~~~~~~~~~~~ 09-03-2021 Read .obj file into memory - vertices and faces ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<unsigned int>> faces;
    read(file_in,vertices,faces);
  //~~~~~~~~~~~~~~~~~~~~~ 09-03-2021 Read .obj file into memory - vertices and faces ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

std::map<unsigned int, std::vector<double>> v_index;

for (unsigned int i = 0; i < vertices.size(); i++){
    v_index.emplace(i,vertices[i]);
}
/* WORKS

    */
    printDCEL(D);

    D.vertices();

  // 2. group the triangles into meshes,
  
  // 3. determine the correct orientation for each mesh and ensure all its triangles 
  //    are consistent with this correct orientation (ie. all the triangle normals 
  //    are pointing outwards).
  
  // 4. merge adjacent triangles that are co-planar into larger polygonal faces.


  // 5. write the meshes with their faces to a valid CityJSON output file. ~~~~~~~~~~~~~~ 09-03-2021~~~~~~~~~~~~~//
    std::fstream fl;
    fl.open (file_out,std::fstream::in | std::fstream::out | std::fstream::trunc);
    fl<<"{\n\"type\": \"CityJSON\",\n\"version\": \"1.0\",\n";
    fl << "\"CityObjects\": {\"id-1\" : {\n\t\"type\": \"Building\",\n\t\"geometry\": [{\n\t\t\"type\": \"MultiSurface\",\n\t\t\"lod\": 2,\n\t\t\"boundaries\": [\n\t\t\t";



    for (auto const& i : D.faces()) {
        unsigned int origin = i->exteriorEdge->origin->i; //->exteriorEdge->origin->i;
        unsigned int destination = i->exteriorEdge->destination->i;
        unsigned int previous = i->exteriorEdge->prev->origin->i;
        if (i == D.faces().back()) {fl << "[[" << origin <<", "<<destination<<", "<<previous << "]]\n\t\t]\n\t}]\n}},\n"; break;}
        fl << "[[" << origin <<", "<<destination<<", "<<previous << "]], ";
    }
    fl << "\"vertices\": [\n";

    for (auto const& i : vertices) {
        double x = i[0]; //->exteriorEdge->origin->i;
        double y = i[1];
        double z = i[2];
        if (i == vertices.back()) {fl << "\t[" << x <<", "<< y <<", "<< z << "]\n\t]\n}"; break;}
        fl << "\t[" << x <<", "<< y<<", "<< z << "],\n";
    }

    fl.close();
  // 5. write the meshes with their faces to a valid CityJSON output file. ~~~~~~~~~~~~~~ 09-03-2021~~~~~~~~~~~~~//

  return 0;
}


void printDCEL(DCEL & D) {

  // Quick check if there is an invalid element
  auto element = D.findInValid();
  if ( element == nullptr ) {
    // Beware that a 'valid' DCEL here only means there are no dangling links and no elimated elements.
    // There could still be problems like links that point to the wrong element.
    std::cout << "DCEL is valid\n";
  } else {
    std::cout << "DCEL is NOT valid ---> ";
    std::cout << *element << "\n";
  }

  // iterate all elements of the DCEL and print the info for each element
  const auto & vertices = D.vertices();
  const auto & halfEdges = D.halfEdges();
  const auto & faces = D.faces();
  std::cout << "DCEL has:\n";
  std::cout << " " << vertices.size() << " vertices:\n";
  for ( const auto & v : vertices ) {
    std::cout << "  * " << *v << "\n";
  }
  std::cout << " " << halfEdges.size() << " half-edges:\n";
  for ( const auto & e : halfEdges ) {
    std::cout << "  * " << *e << "\n";
  }
  std::cout << " " << faces.size() << " faces:\n";
  for ( const auto & f : faces ) {
    std::cout << "  * " << *f << "\n";
  }

}


void DemoDCEL() {

  std::cout << "/// STEP 1 Creating empty DCEL...\n";
  DCEL D;
  printDCEL(D);

  /*

  v2 (0,1,0)
   o
   |\
   | \
   |  \
   o---o v1 (1,0,0)
  v0
  (0,0,0)

  We will construct the DCEL of a single triangle
  in the plane z=0 (as shown above).

  This will require:
    3 vertices
    6 halfedges (2 for each edge)
    1 face

  */
  std::cout << "\n/// STEP 2 Adding triangle vertices...\n";
  Vertex* v0 = D.createVertex(0,0,0,1);
  Vertex* v1 = D.createVertex(1,0,0,2);
  Vertex* v2 = D.createVertex(0,1,0,3);
  printDCEL(D);

  std::cout << "\n/// STEP 3 Adding triangle half-edges...\n";
  HalfEdge* e0 = D.createHalfEdge();
  HalfEdge* e1 = D.createHalfEdge();
  HalfEdge* e2 = D.createHalfEdge();
  HalfEdge* e3 = D.createHalfEdge();
  HalfEdge* e4 = D.createHalfEdge();
  HalfEdge* e5 = D.createHalfEdge();
  printDCEL(D);

  std::cout << "\n/// STEP 4 Adding triangle face...\n";
  Face* f0 = D.createFace();
  printDCEL(D);

  std::cout << "\n/// STEP 5 Setting links...\n";
  e0->origin = v0;
  e0->destination = v1;
  e0->twin = e3;
  e0->next = e1;
  e0->prev = e2;
  e0->incidentFace = f0;

  e3->origin = v1;
  e3->destination = v0;
  e3->twin = e0;
  e3->next = e5;
  e3->prev = e4;

  /*
  If a half-edge is incident to 'open space' (ie not an actual face with an exterior boundary),
  we use the infiniteFace which is predifined in the DCEL class
  */
  e3->incidentFace = D.infiniteFace();

  e1->origin = v1;
  e1->destination = v2;
  e1->twin = e4;
  e1->next = e2;
  e1->prev = e0;
  e1->incidentFace = f0;

  e4->origin = v2;
  e4->destination = v1;
  e4->twin = e1;
  e4->next = e3;
  e4->prev = e5;
  e4->incidentFace = D.infiniteFace();

  e2->origin = v2;
  e2->destination = v0;
  e2->twin = e5;
  e2->next = e0;
  e2->prev = e1;
  e2->incidentFace = f0;

  e5->origin = v0;
  e5->destination = v2;
  e5->twin = e2;
  e5->next = e4;
  e5->prev = e3;
  e5->incidentFace = D.infiniteFace();

  f0->exteriorEdge = e0;
  printDCEL(D);


  std::cout << "\n/// STEP 6 Traversing exterior vertices of f0...\n";
  /*
  if all is well in the DCEL, following a chain of half-edges (ie keep going to e.next)
  should lead us back the the half-edge where we started.
  */
  HalfEdge* e = f0->exteriorEdge;
  const HalfEdge* e_start = e;
  do {
    std::cout << " -> " << *e->origin << "\n";
    e = e->next;
  } while ( e_start!=e) ; // we stop the loop when e_start==e (ie. we are back where we started)


  std::cout << "\n/// STEP 7 eliminating v0...\n";
  v0->eliminate();
  printDCEL(D);

  /*
  We just eliminated v0. At the same time we know there are elements that still
  pointers to v0 (ie the edges e0, e2, e3, e5). This means we can NOT call D.cleanup()!
  If you do this anyways, the program may crash.

  Eg. if you uncomment the following there could be a crash/stall of the program.
  */
  // D.cleanup(); // this will remove v0 from memory (because we just eliminated v0 and the cleanup() function simply removes all the eliminated elements)
  // std::cout << *v0 << "\n"; // we try to access that memory, but v0 is gone -> undefined behaviour
  // std::cout << *e0->origin << "\n"; // this equivalent to the previous line (both point to the same memory address)


  std::cout << "\n/// STEP 8 eliminating all the remaining DCEL elements\n";
  for ( const auto & v : D.vertices() ) {
    v->eliminate();
  }
  for ( const auto & e : D.halfEdges() ) {
    e->eliminate();
  }
  for ( const auto & f : D.faces() ) {
    f->eliminate();
  }
  printDCEL(D);

  std::cout << "\n/// STEP 9 cleaning up the DCEL\n";
  D.cleanup();
  printDCEL(D);

}

#pragma clang diagnostic pop