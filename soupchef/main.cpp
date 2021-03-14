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
#include <deque>
#include <algorithm>


#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmultichar"
// forward declarations; these functions are given below main()
void DemoDCEL();
void printDCEL(DCEL & D);
/*bool rayTriangleIntersect(const std::vector<double> orig, std::vector<double> dir, const Face* f)
{
    auto v0= f->exteriorEdge->origin;
    auto v1= f->exteriorEdge->destination;
    auto v2= f->exteriorEdge->prev->origin;
    std::vector<double> v0v1 = {v1->x - v0->x, v1->y-v0->y, v1->z-v0->z};
    std::vector<double> v0v2 = {v2->x - v0->x, v2->y-v0->y, v2->z-v0->z};


    std::vector<double> pvec = {dir[1] * v0v2[2] - dir[2]*v0v2[1], -(dir[0]*v0v2[2]-dir[2]*v0v2[0]), dir[0]*v0v2[1-dir[1]*v0v2[0]]};
            //dir.cross(v0v2);
    float det = v0v1[0]*pvec[0] + v0v1[1]*pvec[1] + v0v1[2]*pvec[2];
            //v0v1.dot(pvec);

    if (det < 0.000001)
        return false;

    float invDet = 1.0 / det;

    std::vector<double> tvec = {orig[0] - v0->x, orig[1]-v0->y, orig[2]-v0->z};

    float tvec_dot_pvec = tvec[0]*pvec[0] + tvec[1]*pvec[1] + tvec[2] * pvec[2];
    float u = tvec_dot_pvec * invDet;

    if (u < 0 || u > 1)
        return false;

    std::vector<double> qvec =  {tvec[1] * v0v1[2] - tvec[2]*v0v1[1], -(tvec[0]*v0v1[2]-tvec[2]*v0v1[0]), tvec[0]*v0v1[1-tvec[1]*v0v1[0]]};
            //tvec.cross(v0v1);

    float dir_dot_qvec = dir[0]*qvec[0] + dir[1]*qvec[1] + dir[2] * qvec[2];
    float v = dir_dot_qvec * invDet;

    if (v < 0 || u + v > 1)
        return false;

    float v0v2_dot_qvec = v0v2[0]*qvec[0] + v0v2[1]*qvec[1] + v0v2[2] * qvec[2];
    float t = v0v2_dot_qvec * invDet;

    if (t > 0.000001) // ray intersection
    {
        return true;
    }
    else // This means that there is a line intersection but not a ray intersection.
        return true;

    return true;
}
*/ // Other intersection method. Not used.
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
    std::cout<<dot<<" "<<angle <<std::endl;
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


    std::cout<<"flipped\n";
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
        o = {minc[0] - 5, ((maxc[1] - minc[1]) / 2) + minc[1], ((maxc[2] - minc[2]) / 1.5) + minc[2]};
        d = {mesh_faces.back()->exteriorEdge->origin->x,
             mesh_faces.back()->exteriorEdge->origin->y,
             mesh_faces.back()->exteriorEdge->origin->z}; //destination of ray -> a vertex of the last face of the mesh.


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
void mergeCoPlanarFaces(DCEL & D) {
    for (auto const &mesh : D.infiniteFace()->holes) {
        auto mesh_faces = getFaces(mesh);

    std::vector<Face* > traversed_faced;
    for (auto const& f : mesh_faces){
        if (f->isEliminated()){continue;}
        traversed_faced.push_back(f);
        std::vector<double> curr_norm = Normal(f);
        std::vector<HalfEdge*> ne012 = {f->exteriorEdge->twin, f->exteriorEdge->next->twin, f->exteriorEdge->prev->twin};
        for (auto const& neigh : ne012) {
            if (std::find(traversed_faced.begin(), traversed_faced.end(), neigh->incidentFace) !=
                traversed_faced.end()) { continue; }
            if (neigh->isEliminated()){continue;}
            std::vector<double> neigh_norm = Normal(neigh->incidentFace);
            if (angle(curr_norm,neigh_norm) == 0) {
                neigh->twin->twin->eliminate(); // Eliminate the neighboring edge
                neigh->twin->eliminate(); //Eliminate the twin of the neigboring edge. aka the edge of current f.
                neigh->incidentFace->eliminate();
                traversed_faced.push_back(neigh->incidentFace);

                HalfEdge *e = neigh->incidentFace->exteriorEdge;
                const HalfEdge *e_start = e;
                do {
                    e->incidentFace = neigh->twin->incidentFace;
                    e = e->next;
                } while (e_start != e); // we stop the loop when e_start==e (ie. we are back where we started)

               // neigh->next->incidentFace = neigh->twin->incidentFace; // Set the incident face of the next edge of the coplanar neighbor as the current face
               // neigh->prev->incidentFace = neigh->twin->incidentFace;// Set the incident face of the previous edge of the coplanar neighbor as the current face
                if (f->hasDanglingLink()){f->exteriorEdge = f->exteriorEdge->prev;}
                neigh->incidentFace->exteriorEdge = D.halfEdges().front().get();

                neigh->twin->next->prev = neigh->prev;
                neigh->twin->prev->next = neigh->next;

                neigh->prev->next = neigh->twin->next;
                neigh->next->prev = neigh->twin->prev;

                neigh->twin = f->exteriorEdge->twin;
                neigh->twin->twin = f->exteriorEdge;
                int a=0;
            }
        }

    }
    printDCEL(D);

    for (auto const& edge : D.halfEdges()){
    if (edge->hasDanglingLink()){
        if (edge->twin->isEliminated()){edge.get()->twin = edge->next;}
        if (edge->next->isEliminated()){edge.get()->next = edge->prev;}
        if (edge->prev->isEliminated()){edge.get()->prev = edge->twin;}
    }
    if (edge->incidentFace->isEliminated()){
        edge.get()->incidentFace = edge->next->incidentFace;
    }
}
        for (auto const& face : D.faces()) {
            if (face->hasDanglingLink()) {
                if (face->exteriorEdge->isEliminated()){face->exteriorEdge=face->exteriorEdge->next;}
                else if (face->exteriorEdge->next->isEliminated()){face->exteriorEdge->next = face->exteriorEdge->prev;}
                else if (face->exteriorEdge->prev->isEliminated()){face->exteriorEdge->prev = face->exteriorEdge;}
            }
        }
    printDCEL(D);

//2nd layer
        for (auto const& edge : D.halfEdges()){
            if (edge->hasDanglingLink()){
                edge->twin=D.halfEdges().front().get();
            int ad=0;
            }    }


    //D.cleanup();
        printDCEL(D);
    std::map<Face*, std::vector<HalfEdge*>> FtoE;
    for (auto const& face : D.faces()){FtoE[face.get()];}
    for (auto const& edge : D.halfEdges()){
        FtoE[edge->incidentFace].push_back(edge.get());
    }

}}
// 5.
void exportCityJSON(DCEL & D, const char *file_out) {
  // to do
}


int main(int argc, const char * argv[]) {
    const char *file_in = "/home/konstantinos/Desktop/TUDelft-Courses/Q3/GEO1004/hw2/bk_soup.obj";
    const char *file_out = "/home/konstantinos/Desktop/TUDelft-Courses/Q3/GEO1004/hw2/bk.json";




    // Demonstrate how to use the DCEL to get you started (see function implementation below)
    // you can remove this from the final code
    //DemoDCEL();

    // create an empty DCEL
    DCEL D;
    // 1. read the triangle soup from the OBJ input file and convert it to the DCEL,
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<unsigned int>> faces;
    std::unordered_map<std::pair<unsigned int,unsigned int>, HalfEdge *, boost::hash<std::pair<unsigned int, unsigned int>>> hashmap2;
    importOBJ(D, file_in, vertices, faces, hashmap2);
    printDCEL(D);

    // 2. group the triangles into meshes,
    auto mesh_vertices = groupTriangles(D);
    // 3. determine the correct orientation for each mesh and ensure all its triangles
    //    are consistent with this correct orientation (ie. all the triangle normals
    //    are pointing outwards).

    orientMeshes(D, mesh_vertices);
    printDCEL(D);
    // 4. merge adjacent triangles that are co-planar into larger polygonal faces.
    mergeCoPlanarFaces(D);
    printDCEL(D);
    D.cleanup();
    printDCEL(D);
    //printDCEL(D); 5. write the meshes with their faces to a valid CityJSON output file. ~~~~~~~~~~~~~~ 09-03-2021~~~~~~~~~~~~~//
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
            fl << "[[";
            HalfEdge *e = i->exteriorEdge;
            const HalfEdge *e_start = e;
            do {
                if (e->next != e_start) {
                fl << e->origin->i-1 << ",";
                } else {fl << e->origin->i-1;}
                e = e->next;
            } while (e_start != e); // we stop the loop when e_start==e (ie. we are back where we started)
            fl << "]]\n\t\t]\n\t}]\n}},\n";
            break;
        }
        fl << "[[";
        HalfEdge *e = i->exteriorEdge;
        const HalfEdge *e_start = e;
        do {
            if (e->next != e_start) {
                fl << e->origin->i-1 << ",";
            } else {fl << e->origin->i-1;}

            e = e->next;
        } while (e_start != e); // we stop the loop when e_start==e (ie. we are back where we started)
        fl << "]],";
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

    //DemoDCEL();
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