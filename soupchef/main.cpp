#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include<map>
#include <unordered_map>
#include "DCEL.hpp"
#include <cmath>
#include <boost/functional/hash.hpp>

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
std::vector<double> cornerpoints(std::vector<std::vector<double>> v, const std::string& minmax){
    if (minmax == "max"){
        float maxx = -2000.0;
        float maxy = -2000.0;
        float maxz = -2000.0;
        for (auto & i : v) {
            if (i[0] >= maxx) maxx = i[0];
            if (i[1] >= maxy) maxy = i[1];
            if (i[2] >= maxz) maxz = i[2];
        }
        return {maxx+1,maxy+1,maxz+1}; //+1 to go out of bounds
    }
    else if (minmax == "min") {
        float minx = 2000.0;
        float miny = 2000.0;
        float minz = 2000.0;
        for (auto & i : v) {
            if (i[0] <= minx) minx = i[0];
            if (i[1] <= miny) miny = i[1];
            if (i[2] <= minz) minz = i[2];
        }
        return {minx-1,miny-1,minz-1}; //-1 to go out of bounds
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~  19/02/21 Function to calculate SIGNED VOLUME of tetrahedron  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
float signed_volume(Vertex* a, Vertex* b, Vertex* c, std::vector<double> d) {

    std::vector<double> ad = {a->x-d[0], a->y-d[1], a->z-d[2]};
    std::vector<double> bd = {b->x-d[0], b->y-d[1], b->z-d[2]};
    std::vector<double> cd = {c->x-d[0], c->y-d[1], c->z-d[2]};
    std::vector<double> crossbdcd = {bd[1]*cd[2]-bd[2]*cd[1], -(bd[0]*cd[2]-bd[2]*cd[0]), bd[0]*cd[1]-bd[1]*cd[0]};
    auto dotadcross = ad[0]*crossbdcd[0] + ad[1]*crossbdcd[1] + ad[2]*crossbdcd[2];
    return dotadcross/6;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~  19/02/21 Function to calculate SIGNED VOLUME of tetrahedron  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//~~~~~~~~~~~~~~~~~~~~~~~~~~~  19/02/21 Function to fine plane equation and test if 2 two points are on the opposite side ot the plane  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
bool planetest(std::vector<double> orig, std::vector<double> des, Vertex* v, Vertex* vop1, Vertex* vop2 ){
    double a = (des[1] - orig[1]) * (v->z - orig[2]) - (v->y - orig[1]) * (des[2] - orig[2]);
    double b = (v->x - orig[0]) * (des[2] - orig[2]) - (des[0] - orig[0]) * (v->z - orig[2]);
    double c = (des[0] - orig[0]) * (v->y - orig[1]) - (des[1] - orig[1]) * (v->x - orig[0]);
    double d = (-a * orig[0] - b * orig[1] - c * orig[2]);
    double dist1 = ((a * vop1->x) + (b * vop1->y) + (c * vop1->z) + d) / sqrt((a * a) + (b * b) + (c * c));
    double dist2 = ((a * vop2->x) + (b * vop2->y) + (c * vop2->z) + d) / sqrt((a * a) + (b * b) + (c * c));
    if (dist1 * dist2 > 0) {return false;} //~~~~~~~~~~21/02/2021 ">="disregard the case where segment lies on the triangle plane ~~~~~ NOTE: fixes major issues but ignores some valid cases ~~~~~~~~//
    return true;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~  19/02/21 Function to fine plane equation and test if 2 two points are on the opposite side ot the plane  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//~~~~~~~~~~~~~~~~~~~~~~~~~~~  19/02/21 Function to test intersection between line segment and triangle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
bool intersects(std::vector<double> o, std::vector<double> d, Face* i) {

    float signed_or = signed_volume(i->exteriorEdge->origin, i->exteriorEdge->destination, i->exteriorEdge->prev->origin, o);
    float signed_dest = signed_volume(i->exteriorEdge->origin, i->exteriorEdge->destination, i->exteriorEdge->prev->origin, d);
    if (signed_or * signed_dest > 0) {
        return false;
    } else {
        bool inter_v0 = planetest(o, d, i->exteriorEdge->origin, i->exteriorEdge->destination, i->exteriorEdge->prev->origin);
        bool inter_v1 = planetest(o, d, i->exteriorEdge->destination, i->exteriorEdge->origin, i->exteriorEdge->prev->origin);
        bool inter_v2 = planetest(o, d, i->exteriorEdge->prev->origin, i->exteriorEdge->origin, i->exteriorEdge->destination);
        if (inter_v0 && inter_v1 && inter_v2) { return true; }
        return false;
    }
}

Face* min_distance(std::vector<double> o, std::vector<Face *> ray_face){
    std::map<Face*,double> dist;
    for(auto const &j: ray_face){
        std::vector<double> pointx3dist;
        double d0 = sqrt(pow(o[0]-j->exteriorEdge->origin->x,2.0)+pow(o[1]-j->exteriorEdge->origin->y,2.0)+ pow(o[2]-j->exteriorEdge->origin->z,2.0));
        pointx3dist.push_back(d0);
        double d1 = sqrt(pow(o[0]-j->exteriorEdge->destination->x,2.0)+pow(o[1]-j->exteriorEdge->destination->y,2.0)+ pow(o[2]-j->exteriorEdge->destination->z,2.0));
        pointx3dist.push_back(d1);
        double d2 = sqrt(pow(o[0]-j->exteriorEdge->prev->destination->x,2.0)+pow(o[1]-j->exteriorEdge->prev->destination->y,2.0)+ pow(o[2]-j->exteriorEdge->prev->destination->z,2.0));
        pointx3dist.push_back(d2);
        double min = d0;

        for(unsigned int i=0; i<pointx3dist.size(); i++){
            if (pointx3dist[i]<=min){min=pointx3dist[i];}
        }
        dist.insert(std::pair(j,min));
    }
    auto min = dist.begin();
    Face * in;
    for(auto const &k: dist){
        if (k.second<=min->second){
            in = k.first;
        }
    }
return in;
}

std::vector<double> Normal(std::vector<double> v0,std::vector<double> v1, std::vector<double> v2){
    std::vector<double> U = {v1[0]-v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    std::vector<double> V = {v2[0]-v0[0], v2[1] - v0[1], v2[2] - v0[2]};
    std::vector<double> N = {U[1]*V[2] - U[2]*V[1],
                             U[2]*V[0] - U[0]*V[2],
                             U[0]*V[1] - U[1]*V[0]};
    double mag = sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
    const std::vector<double> Nor = {N[0]/mag, N[1]/mag, N[2]/mag};
    return(Nor);
}

/* 
  Example functions that you could implement. But you are 
  free to organise/modify the code however you want.
  After each function you should have a DCEL without invalid elements!
*/
// 1.
void importOBJ(DCEL & D, const char *file_in,std::vector<std::vector<double>>& vertices,std::vector<std::vector<unsigned int>>& faces) {

    std::unordered_map<unsigned int, Vertex *> hashmap;
    read(file_in, vertices, faces);
    unsigned int c = 0;
    for (auto const &i : vertices) {
        c++;
        Vertex *v = D.createVertex(i[0], i[1], i[2], c);
        hashmap[c] = v;
    }
    std::cout<<hashmap[8] <<std::endl;
    std::unordered_map<std::pair<unsigned int,unsigned int>, HalfEdge *, boost::hash<std::pair<unsigned int, unsigned int>>> hashmap2;

    for (auto const &j: faces) {
        Face *f = D.createFace();
        HalfEdge *e0 = D.createHalfEdge();
        HalfEdge *e1 = D.createHalfEdge();
        HalfEdge *e2 = D.createHalfEdge();

        e0->origin = hashmap[j[0]]; //j[0] = origin vertex index
        e0->destination = hashmap[j[1]]; //j[0] = destination vertex index
        e0->next = e1;
        e0->prev = e2;
        e0->incidentFace = f;
        std::pair<unsigned int, unsigned int> p;
        if (j[0]<=j[1]){p.first = j[0]; p.second = j[1];}
        else {p.first = j[1]; p.second = j[0];}

        if (hashmap2.count(p)==0) hashmap2[p] = e0;
        else {
            HalfEdge * etemp = hashmap2[p];
            e0->twin = etemp;
            hashmap2[p]->twin = e0;
        }

        e1->origin = hashmap[j[1]]; //j[0] = origin vertex index
        e1->destination = hashmap[j[2]]; //j[0] = destination vertex index
        e1->next = e2;
        e1->prev = e0;
        e1->incidentFace = f;
        std::pair<unsigned int, unsigned int> p1;
        if (j[1]<=j[2]){p1.first = j[1]; p1.second = j[2];}
        else {p1.first = j[2]; p1.second = j[1];}

        if (hashmap2.count(p1)==0) hashmap2[p1] = e1;
        else {
            HalfEdge * etemp1 = hashmap2[p1];
            e1->twin = etemp1;
            hashmap2[p1]->twin = e1;
        }

        e2->origin = hashmap[j[2]]; //j[0] = origin vertex index
        e2->destination = hashmap[j[0]]; //j[0] = destination vertex index
        e2->next = e0;
        e2->prev = e1;
        e2->incidentFace = f;
        std::pair<unsigned int, unsigned int> p2;
        if (j[2]<=j[0]){p2.first = j[2]; p2.second = j[0];}
        else {p2.first = j[0]; p2.second = j[2];}

        if (hashmap2.count(p2)==0) hashmap2[p2] = e2;
        else {
            HalfEdge * etemp2 = hashmap2[p2];
            e2->twin = etemp2;
            hashmap2[p2]->twin = e2;
        }

        f->exteriorEdge = e0;
        D.infiniteFace()->holes.push_back(e0);





    }
}
// 2.
void groupTriangles(DCEL & D) {
  // to do
}
// 3.
void orientMeshes(DCEL & D, std::vector<double> o, std::vector<double> d ,std::vector<std::vector<double>>& vertice) {
  o = cornerpoints(vertice,"min"); //origin of ray
  d = cornerpoints(vertice,"max"); //destination of ray
  std::vector<Face *> ray_face;
  for(auto const& i : D.faces()){
      if (intersects(o,d,i.get())){
          ray_face.push_back(i.get());
      }
  }
  Face* nearest = min_distance(o, ray_face);

  std::vector<double> fv0 {nearest->exteriorEdge->origin->x, nearest->exteriorEdge->origin->y,nearest->exteriorEdge->origin->z};
  std::vector<double> fv1 {nearest->exteriorEdge->destination->x, nearest->exteriorEdge->destination->y,nearest->exteriorEdge->destination->z};
  std::vector<double> fv2 {nearest->exteriorEdge->prev->origin->x, nearest->exteriorEdge->prev->origin->y,nearest->exteriorEdge->prev->origin->z};
  std::vector<double> e;

  e = Normal(fv0,fv2,fv1);

  //double angle = acos((e[0]*o[0] + e[1]*o[1] + e[2]*o[2]) / sqrt((e[0]*e[0]+e[1]*e[1]+e[2]*e[2]) * (o[0]*o[0]+o[1]*o[1]+o[2]*o[2])));
  double dot = e[0]*o[0] + e[1]*o[1] + e[2]*o[2]; //The correct orientation is when dot is positive
  std::cout<<dot<<std::endl;
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
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<unsigned int>> faces;
    importOBJ(D, file_in, vertices, faces);
    printDCEL(D);

    D.vertices();

    // 2. group the triangles into meshes,

    // 3. determine the correct orientation for each mesh and ensure all its triangles
    //    are consistent with this correct orientation (ie. all the triangle normals
    //    are pointing outwards).
    std::vector<double> ray_origin;
    std::vector<double> ray_destination;
orientMeshes(D,ray_origin,ray_destination, vertices);
    // 4. merge adjacent triangles that are co-planar into larger polygonal faces.


    // 5. write the meshes with their faces to a valid CityJSON output file. ~~~~~~~~~~~~~~ 09-03-2021~~~~~~~~~~~~~//
    std::fstream fl;
    fl.open(file_out, std::fstream::in | std::fstream::out | std::fstream::trunc);
    fl << "{\n\"type\": \"CityJSON\",\n\"version\": \"1.0\",\n";
    fl
            << "\"CityObjects\": {\"id-1\" : {\n\t\"type\": \"Building\",\n\t\"geometry\": [{\n\t\t\"type\": \"MultiSurface\",\n\t\t\"lod\": 2,\n\t\t\"boundaries\": [\n\t\t\t";

//CHECK THIS, IT WORKS WITH THE OLDER IMPLEMENTATION OF INDEXING
    for (auto const &i : D.faces()) {
        unsigned int origin = i->exteriorEdge->origin->i; //->exteriorEdge->origin->i;
        unsigned int destination = i->exteriorEdge->destination->i;
        unsigned int previous = i->exteriorEdge->prev->origin->i;
        if (i == D.faces().back()) {
            fl << "[[" << origin-1 << ", " << destination-1 << ", " << previous-1 << "]]\n\t\t]\n\t}]\n}},\n";
            break;
        }
        fl << "[[" << origin-1 << ", " << destination-1 << ", " << previous-1 << "]], ";
    }
    fl << "\"vertices\": [\n";

    for (auto const &i : vertices) {
        double x = i[0]; //->exteriorEdge->origin->i;
        double y = i[1];
        double z = i[2];
        if (i == vertices.back()) {
            fl << "\t[" << x << ", " << y << ", " << z << "]\n\t]\n}";
            break;
        }
        fl << "\t[" << x << ", " << y << ", " << z << "],\n";
    }

    fl.close();
    // 5. write the meshes with their faces to a valid CityJSON output file. ~~~~~~~~~~~~~~ 09-03-2021~~~~~~~~~~~~~//

    DemoDCEL();
    return 0;
}

    void printDCEL(DCEL &D) {

        // Quick check if there is an invalid element
        auto element = D.findInValid();
        if (element == nullptr) {
            // Beware that a 'valid' DCEL here only means there are no dangling links and no elimated elements.
            // There could still be problems like links that point to the wrong element.
            std::cout << "DCEL is valid\n";
        } else {
            std::cout << "DCEL is NOT valid ---> ";
            std::cout << *element << "\n";
        }

        // iterate all elements of the DCEL and print the info for each element
        const auto &vertices = D.vertices();
        const auto &halfEdges = D.halfEdges();
        const auto &faces = D.faces();
        std::cout << "DCEL has:\n";
        std::cout << " " << vertices.size() << " vertices:\n";
        for (const auto &v : vertices) {
            std::cout << "  * " << *v << "\n";
        }
        std::cout << " " << halfEdges.size() << " half-edges:\n";
        for (const auto &e : halfEdges) {
            std::cout << "  * " << *e << "\n";
        }
        std::cout << " " << faces.size() << " faces:\n";
        for (const auto &f : faces) {
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
        Vertex *v0 = D.createVertex(0, 0, 0, 1);
        Vertex *v1 = D.createVertex(1, 0, 0, 2);
        Vertex *v2 = D.createVertex(0, 1, 0, 3);
        printDCEL(D);

        std::cout << "\n/// STEP 3 Adding triangle half-edges...\n";
        HalfEdge *e0 = D.createHalfEdge();
        HalfEdge *e1 = D.createHalfEdge();
        HalfEdge *e2 = D.createHalfEdge();
        HalfEdge *e3 = D.createHalfEdge();
        HalfEdge *e4 = D.createHalfEdge();
        HalfEdge *e5 = D.createHalfEdge();
        printDCEL(D);

        std::cout << "\n/// STEP 4 Adding triangle face...\n";
        Face *f0 = D.createFace();
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
        HalfEdge *e = f0->exteriorEdge;
        const HalfEdge *e_start = e;
        do {
            std::cout << " -> " << *e->origin << "\n";
            e = e->next;
        } while (e_start != e); // we stop the loop when e_start==e (ie. we are back where we started)


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
        for (const auto &v : D.vertices()) {
            v->eliminate();
        }
        for (const auto &e : D.halfEdges()) {
            e->eliminate();
        }
        for (const auto &f : D.faces()) {
            f->eliminate();
        }
        printDCEL(D);

        std::cout << "\n/// STEP 9 cleaning up the DCEL\n";
        D.cleanup();
        printDCEL(D);

    }

#pragma clang diagnostic pop