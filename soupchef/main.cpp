#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <unordered_map>
#include "DCEL.hpp"
#include <cmath>
#include <boost/functional/hash.hpp>
#include <stack>
#include <algorithm>




#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmultichar"
// forward declarations; these functions are given below main()
void DemoDCEL();
void printDCEL(DCEL & D);
//~~~~~~~~~~~~~~~~~~~~~ 09-03-2021 Read .obj file into memory - vertices and faces ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
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
        return {maxx,maxy,maxz}; //+1 to go out of bounds
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
        return {minx,miny,minz}; //-1 to go out of bounds
    }
}
float signed_volume(Vertex* a, Vertex* b, Vertex* c, std::vector<double> d) {

    std::vector<double> ad = {a->x-d[0], a->y-d[1], a->z-d[2]};
    std::vector<double> bd = {b->x-d[0], b->y-d[1], b->z-d[2]};
    std::vector<double> cd = {c->x-d[0], c->y-d[1], c->z-d[2]};
    std::vector<double> crossbdcd = {bd[1]*cd[2]-bd[2]*cd[1], -(bd[0]*cd[2]-bd[2]*cd[0]), bd[0]*cd[1]-bd[1]*cd[0]};
    auto dotadcross = ad[0]*crossbdcd[0] + ad[1]*crossbdcd[1] + ad[2]*crossbdcd[2];
    return dotadcross/6;
}
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
bool checkNormal(std::vector<double> v0,std::vector<double> v1, std::vector<double> v2, std::vector<double>& o,std::vector<double>& d){
    std::vector<double> U = {v1[0]-v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    std::vector<double> V = {v2[0]-v0[0], v2[1] - v0[1], v2[2] - v0[2]};
    std::vector<double> N = {U[1]*V[2] - U[2]*V[1],
                             U[2]*V[0] - U[0]*V[2],
                             U[0]*V[1] - U[1]*V[0]};
    double mag = sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
    const std::vector<double> Nor = {N[0]/mag, N[1]/mag, N[2]/mag};
    std::vector<double> ray = {d[0]-o[0], d[1]-o[1], d[2]-o[2]};
    double angle = acos((Nor[0]*ray[0] + Nor[1]*ray[1] + Nor[2]*ray[2]) / sqrt((Nor[0]*Nor[0]+Nor[1]*Nor[1]+Nor[2]*Nor[2]) * (ray[0]*ray[0]+ray[1]*ray[1]+ray[2]*ray[2])));
    double dot = Nor[0]*ray[0] + Nor[1]*ray[1] + Nor[2]*ray[2];
    if (dot<0) return true; //The correct orientation is when dot is negative (means that the vectors have opposite-like direction.
    else if (dot>0) return false;
    else std::cout << "dot is 0";  //Think about the vertical case of dot=0/ Think of many cases if it waorks every time.
}
std::vector<double> Normal(const Face* f){
    std::vector<double>  v0 = {f->exteriorEdge->origin->x ,f->exteriorEdge->origin->y,f->exteriorEdge->origin->z};
    std::vector<double>  v1 = {f->exteriorEdge->destination->x ,f->exteriorEdge->destination->y,f->exteriorEdge->destination->z};
    std::vector<double>  v2 = {f->exteriorEdge->prev->origin->x ,f->exteriorEdge->prev->origin->y,f->exteriorEdge->prev->origin->z};
    std::vector<double> U = {v1[0]-v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    std::vector<double> V = {v2[0]-v0[0], v2[1] - v0[1], v2[2] - v0[2]};
    std::vector<double> N = {U[1]*V[2] - U[2]*V[1],
                             U[2]*V[0] - U[0]*V[2],
                             U[0]*V[1] - U[1]*V[0]};
    double mag = sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
    const std::vector<double> Nor = {N[0]/mag, N[1]/mag, N[2]/mag};
    return Nor;
}
double angle(std::vector<double> n1, std::vector<double> n2) {
    double angle = acos((n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) /
                        sqrt((n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]) *
                             (n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2])));
    return angle = angle*180/M_PI;
}
void flip012_201(Face* triangle){
    HalfEdge* e0 = triangle->exteriorEdge;
    HalfEdge* e1 = e0->next;
    HalfEdge* e2 = e0->prev;
    Vertex* t0 = e0->destination;
    Vertex* t1 = e0->origin;
    Vertex* t2 = e0->next->destination;
    //HalfEdge* twinnext = triangle->exteriorEdge->prev->twin;
    //HalfEdge* twinprev = triangle->exteriorEdge->next->twin;

    triangle->exteriorEdge->origin = t0;
    triangle->exteriorEdge->destination = t1 ;
    triangle->exteriorEdge->next = e2;
    triangle->exteriorEdge->prev = e1;

    triangle->exteriorEdge->next->origin = t1;
    triangle->exteriorEdge->next->destination = t2 ;
    triangle->exteriorEdge->next->next = e1;
    triangle->exteriorEdge->next->prev = e0;

    triangle->exteriorEdge->prev->origin = t2;
    triangle->exteriorEdge->prev->destination = t0 ;
    triangle->exteriorEdge->prev->next = e0;
    triangle->exteriorEdge->prev->prev = e2;
}
std::vector<double> Centroid(HalfEdge *f){
    auto v0 = f->origin;
    auto v1 = f->destination;
    auto v2 = f->prev->origin;
    std::vector<double> c = {(v0->x+v1->x+v2->x)/3, (v0->y+v1->y+v2->y)/3, (v0->z+v1->z+v2->z)/3};
    return c;
}
std::vector<Face* > getFaces(HalfEdge* edge){
    std::vector<Face* > toreturn;
    std::stack<Face* > facestack;
    std::vector<Face* > traversed_faced;
    facestack.push(edge->incidentFace);
    while (!facestack.empty()) {
        auto validface = facestack.top();
        traversed_faced.push_back(validface);
        facestack.pop();
        //if (validface->isEliminated()){continue;} //16-03
        std::vector<HalfEdge *> e012 = {validface->exteriorEdge, validface->exteriorEdge->next,
                                        validface->exteriorEdge->prev};

        for (auto const &edgex : e012) {
            if (std::find(traversed_faced.begin(), traversed_faced.end(), edgex->twin->incidentFace) !=
                traversed_faced.end()) {
                continue;
            }
            facestack.push(edgex->twin->incidentFace);
            traversed_faced.push_back(edgex->twin->incidentFace);
        }
        toreturn.push_back(validface);
    }
    return toreturn;
}
// 1.
void importOBJ(DCEL & D, const char *file_in,std::vector<std::vector<double>>& vertices,std::vector<std::vector<unsigned int>>& faces,std::unordered_map<std::pair<unsigned int,unsigned int>, HalfEdge *, boost::hash<std::pair<unsigned int, unsigned int>>>& hashmap2) {

    std::unordered_map<unsigned int, Vertex *> hashmap;
    read(file_in, vertices, faces);
    unsigned int c = 0;
    for (auto const &i : vertices) {
        c++;
        Vertex *v = D.createVertex(i[0], i[1], i[2], c);
        hashmap[c] = v;
    }

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
    }
}
// 2.
std::vector<std::vector<std::vector<double>>> groupTriangles(DCEL & D) {
    std::vector<Vertex*> verts;
    std::vector<std::vector<double>> v_toreturn;
    std::stack<Face* > facestack;
    std::vector<Face* > traversed_faced;
    std::vector<std::vector<std::vector<double>>> mesh_vertices;

    for (auto const& face : D.faces()) {
        if (std::find(traversed_faced.begin(),traversed_faced.end(), face.get()) != traversed_faced.end())
        {continue;}
        facestack.push(face.get());
        verts={};
        v_toreturn={};

        while (!facestack.empty()) {
            auto validface = facestack.top();
            traversed_faced.push_back(validface);
            facestack.pop();
            std::vector<HalfEdge *> e012 = {validface->exteriorEdge, validface->exteriorEdge->next,
                                            validface->exteriorEdge->prev};

            for (auto const &edgex : e012) {
                if (std::find(traversed_faced.begin(), traversed_faced.end(), edgex->twin->incidentFace) !=
                    traversed_faced.end()) {
                    continue;
                }
                facestack.push(edgex->twin->incidentFace);
                traversed_faced.push_back(edgex->twin->incidentFace);

            }
            std::vector<Vertex *> v012 = {validface->exteriorEdge->origin, validface->exteriorEdge->destination,validface->exteriorEdge->prev->origin};
            for (auto const& vertices : v012)
            if (std::find(verts.begin(),verts.end(), vertices) != verts.end())
            {continue;}
            else {verts.push_back(vertices); v_toreturn.push_back({vertices->x,vertices->y,vertices->z});}
        }

        mesh_vertices.push_back(v_toreturn);
        D.infiniteFace()->holes.push_back(face->exteriorEdge);
    }
    return mesh_vertices;
}
// 3.
void orientMeshes(DCEL & D ,std::vector<std::vector<std::vector<double>>>& vertice) {
    unsigned int mesh_count=0;
    for (auto const &mesh : D.infiniteFace()->holes) {
        std::vector<double> o; //origin of ray
        std::vector<double> d; //destination of ray
        auto mesh_faces = getFaces(mesh); //mesh_faces contain only those faces that difine the specific mesh.

        auto minc = cornerpoints(vertice[mesh_count], "min");
        auto maxc = cornerpoints(vertice[mesh_count], "max");
        o = {-10.0 * abs(round(minc[0])), -10.0 * abs(round(minc[1])), minc[2]+(maxc[2] - minc[2])/1.4};
        d = Centroid(mesh_faces.at(mesh_faces.size()-2)->exteriorEdge);//{mesh_faces.back()->exteriorEdge->origin->x,mesh_faces.back()->exteriorEdge->origin->y,mesh_faces.back()->exteriorEdge->origin->z}; //destination of ray -> a vertex of the last face of the mesh.

        std::vector<Face *> ray_face; //triangles that the ray intersects
        for (auto const & i : mesh_faces) {
            if (intersects(o, d, i)) {
                ray_face.push_back(i);
            }
        }
        Face *nearest = min_distance(o, ray_face); //Closest intersecting triangle to ray.

        std::vector<double> fv0{nearest->exteriorEdge->origin->x, nearest->exteriorEdge->origin->y,
                                nearest->exteriorEdge->origin->z};
        std::vector<double> fv1{nearest->exteriorEdge->destination->x, nearest->exteriorEdge->destination->y,
                                nearest->exteriorEdge->destination->z};
        std::vector<double> fv2{nearest->exteriorEdge->prev->origin->x, nearest->exteriorEdge->prev->origin->y,
                                nearest->exteriorEdge->prev->origin->z};

//Check and fix, if necessary, the orientation of the initial face.
//This face is the closest intersecting face of the ray with origin near the middle of bbox for y,z but outside bounds for x and
//destination a random triangle's centroid (to ensure at least one intersecttion is happening).
        if (!checkNormal(fv0, fv1, fv2, o, d))
            flip012_201(nearest);


// Fix the orientation of all the faces of the mesh.  ###########################################################
// Start with "nearest" triangle which has the correct orientation and set the rest based on it. #############
        std::stack<Face *> facestack;
        std::vector<Face *> traversed_faced;
        facestack.push(nearest);
        while (!facestack.empty()) {
            auto validface = facestack.top();
            traversed_faced.push_back(validface);
            facestack.pop();
            std::vector<HalfEdge *> e012 = {validface->exteriorEdge, validface->exteriorEdge->next,
                                            validface->exteriorEdge->prev};

            for (auto const &edgex : e012) {
                if (std::find(traversed_faced.begin(), traversed_faced.end(), edgex->twin->incidentFace) !=
                    traversed_faced.end()) {
                    continue;
                }
                if (edgex->origin == edgex->twin->origin) {
                    flip012_201(edgex->twin->incidentFace);
                }
                facestack.push(edgex->twin->incidentFace);
                traversed_faced.push_back(edgex->twin->incidentFace);
            }
        }
//###############################################################################################################
    mesh_count++;
    }
}
// 4.
void mergeCoPlanarFaces(DCEL & D, DCEL& tempD) {

    for (auto const &mesh : D.infiniteFace()->holes) {
        auto mesh_faces = getFaces(mesh);
      unsigned int c=0;
    std::vector<Face* > traversed_faced;
    for (auto const& face : mesh_faces){

        if (face->isEliminated()){; continue;}
        traversed_faced.push_back(face);
        std::stack<Face *> facestack;
        facestack.push(face);
        std::vector<double> curr_norm = Normal(face);

        while (!facestack.empty()) {
            auto check_face = facestack.top();
            //traversed_faced.push_back(check_face);
            facestack.pop();

            HalfEdge *e = check_face->exteriorEdge;
            const HalfEdge *e_start = e;
            do {
                std::vector<double> neigh_norm = Normal(e->twin->incidentFace);

                if (e->incidentFace == e->twin->incidentFace){

                    e = e->next;
                    c++;
                    continue;
                    }
                if (angle(curr_norm,neigh_norm) <= 1){
                    e->eliminate(); // Eliminate the edge of checking face
                    e->twin->eliminate(); //Eliminate the twin of the checking edge. aka the edge of neighbor face.

                    e->twin->incidentFace->eliminate(); // Eliminate the neighboring face

                    traversed_faced.push_back(e->twin->incidentFace);
                    facestack.push(check_face);

                    //Link redirections to make the merged face valid (NOTE: not the eliminated ones)
                    e->prev->next = e->twin->next;
                    e->twin->next->prev = e->prev;

                    e->next->prev = e->twin->prev;
                    e->twin->prev->next = e->next;

                    e->twin->next->incidentFace = check_face;
                    e->twin->prev->incidentFace = check_face;

                    //e_start=e->next; //prevents infinite loop // when a merge happens code gets out of the while loop an gets back in again to the updated polygonal (now) face.
                    check_face->exteriorEdge = e->next;
                    break;
                }
                e = e->next;
            } while (e_start != e); // we stop the loop when e_start==e (ie. we are back where we started)
        }}
printDCEL(D);


}



}
// 5.
void exportCityJSON(DCEL & D,std::vector<std::vector<double>> vertices, const char *file_out) {
    //5. write the meshes with their faces to a valid CityJSON output file. ~~~~~~~~~~~~~~ 09-03-2021~~~~~~~~~~~~~//
    unsigned int count = 0;
    std::fstream fl;
    fl.open(file_out, std::fstream::in | std::fstream::out | std::fstream::trunc);
    fl << "{\n\"type\": \"CityJSON\",\n\"version\": \"1.0\",\n";
    fl << "\"CityObjects\": {\n";


    std::map<unsigned int,std::vector<Face*>> faces;
    for (auto const&face : D.faces()){
        faces[face->mesh].push_back(face.get());
    }

    for (unsigned int i = 0; i<faces.size(); i++) {

        count++;
        fl << "\"id-" << std::to_string(count)<<"\""
           << ": {\n\t\"type\": \"Building\",\n\t\"geometry\": [{\n\t\t\"type\": \"MultiSurface\",\n\t\t\"lod\": 2,\n\t\t\"boundaries\": [\n\t\t\t";


        for (auto const &i : faces.at(count)) {
            unsigned int origin = i->exteriorEdge->origin->i; //->exteriorEdge->origin->i;
            unsigned int destination = i->exteriorEdge->destination->i;
            unsigned int previous = i->exteriorEdge->prev->origin->i;
            if (i == faces.at(count).back()) {
                fl << "[[";
                HalfEdge *e = i->exteriorEdge;
                const HalfEdge *e_start = e;
                do {
                    if (e->next != e_start) {
                        fl << e->origin->i - 1 << ",";
                    } else { fl << e->origin->i - 1; }
                    e = e->next;
                } while (e_start != e); // we stop the loop when e_start==e (ie. we are back where we started)
                if (i->holes.size() == 0) {
                    if (count == faces.size()) {
                        fl << "]]\n\t\t]\n\t}]\n";
                    }else {fl << "]]\n\t\t]\n\t}]},\n";}
                }
/*Holes*/   else {
                    fl << "],[";
                    HalfEdge *e = i->holes.front();
                    const HalfEdge *e_start = e;
                    do {
                        if (e->next != e_start) {
                            fl << e->origin->i - 1 << ",";
                        } else { fl << e->origin->i - 1; }
                        e = e->next;
                    } while (e_start != e);
                    fl << "]]\n\t\t]\n\t}]\n";
                }
                break;
            }
            fl << "[[";
            HalfEdge *e = i->exteriorEdge;
            const HalfEdge *e_start = e;
            do {
                if (e->next != e_start) {
                    fl << e->origin->i - 1 << ",";
                } else { fl << e->origin->i - 1; }

                e = e->next;
            } while (e_start != e); // we stop the loop when e_start==e (ie. we are back where we started)
            if (i->holes.size() == 0) {
                fl << "]],";
            }
/*Holes*/   else {
                fl << "],[";
                HalfEdge *e = i->holes.front();
                const HalfEdge *e_start = e;
                do {
                    if (e->next != e_start) {
                        fl << e->origin->i - 1 << ",";
                    } else { fl << e->origin->i - 1; }
                    e = e->next;
                } while (e_start != e);
                fl << "]],";
            }
        }






        // 5. write the meshes with their faces to a valid CityJSON output file. ~~~~~~~~~~~~~~ 09-03-2021~~~~~~~~~~~~~//



    }
    fl << "}},\n\"vertices\": [\n";

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
}
    int main(int argc, const char *argv[]) {
        const char *file_in = "/home/konstantinos/Desktop/TUDelft-Courses/Q3/GEO1004/hw2/bk_soup.obj";
        const char *file_out = "/home/konstantinos/Desktop/TUDelft-Courses/Q3/GEO1004/hw2/bk.json";

        DCEL tempD;
        Face *Ftemp = tempD.createFace();
        Vertex *vtemp = tempD.createVertex(0, 0, 0, 1);
        HalfEdge *etemp = tempD.createHalfEdge();

        tempD.faces().front().get()->exteriorEdge = etemp;
        tempD.faces().front().get()->exteriorEdge->origin = vtemp;
        tempD.faces().front().get()->exteriorEdge->destination = vtemp;
        tempD.faces().front().get()->exteriorEdge->next = etemp;
        tempD.faces().front().get()->exteriorEdge->prev = etemp;
        tempD.faces().front().get()->exteriorEdge->twin = etemp;
        tempD.faces().front().get()->exteriorEdge->incidentFace = Ftemp;
        printDCEL(tempD);


        // create an empty DCEL
        DCEL D;
        // 1. read the triangle soup from the OBJ input file and convert it to the DCEL,
        std::vector<std::vector<double>> vertices;
        std::vector<std::vector<unsigned int>> faces;
        std::unordered_map<std::pair<unsigned int, unsigned int>, HalfEdge *, boost::hash<std::pair<unsigned int, unsigned int>>> hashmap2;
        importOBJ(D, file_in, vertices, faces, hashmap2);
        //printDCEL(D);

        // 2. group the triangles into meshes,
        auto mesh_vertices = groupTriangles(D);
        // 3. determine the correct orientation for each mesh and ensure all its triangles
        //    are consistent with this correct orientation (ie. all the triangle normals
        //    are pointing outwards).
        unsigned int c=0;
        for (auto const& mesh : D.infiniteFace()->holes){
            auto mesh_faces = getFaces(mesh);
            c++;
            for (auto const& face : mesh_faces) {
                face->mesh = c;
            }
        }

        orientMeshes(D, mesh_vertices);

        // 4. merge adjacent triangles that are co-planar into larger polygonal faces.
        mergeCoPlanarFaces(D, tempD);

//        printDCEL(D);
        D.cleanup();
//        printDCEL(D);



        // 4.1 search for holes
        for (auto const &Face : D.faces()) {
            std::vector<HalfEdge *> cords;
            HalfEdge *e = Face->exteriorEdge;
            const HalfEdge *e_start = e;
            do {
                // if (e->incidentFace == e->twin->incidentFace && e_start==e){continue;}
                if (e->incidentFace == e->twin->incidentFace) {
                    cords.push_back(e);
                    Face->holes.push_back(e->next);
                    e->eliminate();
                    e->twin->eliminate();

                    e->prev->next = e->twin->next;
                    e->twin->next->prev = e->prev;

                    e->next->prev = e->twin->prev;
                    e->twin->prev->next = e->next;


                    e_start = e->next; //avoid infinite loop;
                }
                e = e->next;
            } while (e_start != e);
        }


//        printDCEL(D);
        D.cleanup();
//        printDCEL(D);

    for (auto const &faces : D.faces()) {
            double distance = 0;double distance2=0;
            if (faces->holes.size() != 0) {

                for (auto const& hole : faces->holes){
                    HalfEdge *e = faces->exteriorEdge;
                    const HalfEdge *e_start = e;
                    do {
                        distance = distance + sqrt(pow((e->destination->x - e->origin->x),2)+pow((e->destination->y - e->origin->y),2)+pow((e->destination->z - e->origin->z),2));
                        e=e->next;
                    } while (e_start != e);
                    HalfEdge *eh = faces->holes.front();
                    const HalfEdge *eh_start = eh;
                    do {
                        distance2 = distance2 + sqrt(pow((eh->destination->x - eh->origin->x),2)+pow((eh->destination->y - eh->origin->y),2)+pow((eh->destination->z - eh->origin->z),2));
                        eh=eh->next;
                    } while (eh_start != eh);
                }
                HalfEdge *interior = faces->exteriorEdge;
                HalfEdge *exterior = faces->holes.front();
                if (distance < distance2){
                    faces->holes.front() = interior;
                    faces->exteriorEdge = exterior;
                }
            }
            }
    for (auto const &e:D.halfEdges()) {
        if (e->hasDanglingLink()) {
            e->incidentFace = e->twin->incidentFace;
        }
    }


        printDCEL(D);
        // 5. Export CityJson file into file_out
        exportCityJSON(D, vertices, file_out);
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

#pragma clang diagnostic pop