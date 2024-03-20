

// magic define number of houdini in pbd_types.h
#define DISTANCE -264586729
#define BEND 5106433


// #include </opt/hfs19.5/houdini/ocl/include/pbd_types.h>
// #include </opt/hfs19.5/houdini/ocl/include/typedefines.h>
// #include </opt/hfs19.5/houdini/ocl/include/matrix.h>
// #include </opt/hfs19.5/houdini/ocl/include/quaternion.h>
// #include </opt/hfs19.5/houdini/ocl/include/reduce.h>

static void distanceUpdateXPBD(
    float timeinc,
    int idx,
    int ptidx,
    global const int *pts,
    global float *L,
    global float *P,
    global const float *pprev,
    global const float *mass,

    float restlength, // restlength will be distance for distance
    float kstiff,
    float kdampratio
    
) {























    // distance will get the two points from pts and calculate the distance between them
    int pt0 = pts[ptidx];
    int pt1 = pts[ptidx + 1];
    float3 p0 = vload3(pt0, P); // get the position of the first point
    float3 p1 = vload3(pt1, P); // get the position of the second point
    float3 prev0 = vload3(pt0, pprev); // get the previous position of the first point
    float3 prev1 = vload3(pt1, pprev); // get the previous position of the second point
    
    float invmass0 = select(0.0f, 1.0f / mass[pt0], mass[pt0] >= 0.0); // get the inverse mass of the first point
    float invmass1 = select(0.0f, 1.0f / mass[pt1], mass[pt1] >= 0.0); // get the inverse mass of the second point
    float gradC_invM_GradC = invmass0 + invmass1;
    if (gradC_invM_GradC == 0.0f) return;

    float3 direct_n = p1 - p0; // get the direction vector between the two points
    float dist_n = length(direct_n); // get the distance between the two points
    if (dist_n < 1e-6f) return;
    if (kstiff == 0.0f) return;

    float l = L[idx * 3]; // get the Lagrange multiplier(also known as lambda)
    float alpha = 1.0 / kstiff / (timeinc * timeinc); // alpha/h^2 XPBD term
    float C = dist_n - restlength; // |x12|-d Constraint term
    float3 direct_n_normalize = direct_n / dist_n; // normalize the direction vector
    float3 gradC = direct_n_normalize; // x12/|x12| Gradient Constraint

    // if kdampratio is 0, then dsum = 0, gamma = 1
    float dsum = 0.0, gamma = 1.0;
    // float beta = kstiff * kdampratio * timeinc * timeinc;
    // float gamma = alpha * beta / timeinc;
    // float dsum = gamma * ((-dot(gradC, p0 - prev0)) + dot(gradC, p1 - prev1));
    // gamma += 1.0f;
    float dL = (-C - alpha * l - dsum) / (gamma * gradC_invM_GradC + alpha);

    float3 dp = -dL * direct_n_normalize;
    p0 += invmass0 * dp;
    p1 += -invmass1 * dp;
    vstore3(p0, pt0, P);
    vstore3(p1, pt1, P);
    L[idx * 3] += dL;

}



static void bendUpdateXPBD(
    float timeinc,
    int idx,
    int ptidx,
    global const int *pts,
    global float *L,
    global float *P,
    global const float *pprev,
    global const float *mass,

    float restlength, // restlength will be angle for bend
    float kstiff,
    float kdampratio
    
) {

    // bend will get the four points of two triangles from pts and calculate the bend between them
    int pt0 = pts[ptidx];
    int pt1 = pts[ptidx + 1];
    int pt2 = pts[ptidx + 2];
    int pt3 = pts[ptidx + 3];
    float3 p0 = vload3(pt0, P); // get the position of the first point
    float3 p1 = vload3(pt1, P); // get the position of the second point
    float3 p2 = vload3(pt2, P); // get the position of the third point
    float3 p3 = vload3(pt3, P); // get the position of the fourth point
    float invmass0 = select(0.0f, 1.0f / mass[pt0], mass[pt0] >= 0.0); // get the inverse mass of the first point
    float invmass1 = select(0.0f, 1.0f / mass[pt1], mass[pt1] >= 0.0); // get the inverse mass of the second point
    float invmass2 = select(0.0f, 1.0f / mass[pt2], mass[pt2] >= 0.0); // get the inverse mass of the third point
    float invmass3 = select(0.0f, 1.0f / mass[pt3], mass[pt3] >= 0.0); // get the inverse mass of the fourth point
    float l = L[idx * 3];

    //     /\ p0
    //    /  \
    // p2/____\p3
    //   \    /
    //    \  /
    //     \/ p1
    float3 n1 = cross(p3 - p0, p2 - p0);
    float3 n2 = cross(p2 - p1, p3 - p1);
    if (length(n1) < 1e-6f || length(n2) < 1e-6f) return;

    float3 n1_normalize_weird = n1 / (length(n1) * length(n1)); // why not use length? (why use square)
    float3 n2_normalize_weird = n2 / (length(n2) * length(n2)); // why not use length? (why use square)
    float3 n1_normalize = normalize(n1);
    float3 n2_normalize = normalize(n2);

    float3 mid_edge = p3 - p2;
    float mid_edge_length = length(mid_edge);
    float inv_mid_edge_length = 1.0f / length(mid_edge);

    float s = select(-1.0f, 1.0f, dot(cross(n1_normalize, n2_normalize), mid_edge) > 0.0);

    float3 gradC0 = s * mid_edge_length * n1_normalize;
    float3 gradC1 = s * mid_edge_length * n2_normalize;
    float3 gradC2 = s * (dot(p0 - p3, mid_edge) * inv_mid_edge_length * n1_normalize 
                        + dot(p1 - p3, mid_edge) * inv_mid_edge_length * n2_normalize);
    float3 gradC3 = s * (dot(p2 - p0, mid_edge) * inv_mid_edge_length * n1_normalize
                        + dot(p2 - p1, mid_edge) * inv_mid_edge_length * n2_normalize);

    float gradC_invM_GradC = invmass0 * dot(gradC0, gradC0) +
                            invmass1 * dot(gradC1, gradC1) +
                            invmass2 * dot(gradC2, gradC2) +
                            invmass3 * dot(gradC3, gradC3);
    if (gradC_invM_GradC == 0.0) return;

    float alpha = 1.0f / kstiff / (timeinc * timeinc);
    float phi = acos(clamp(dot(n1, n2), -1.0f, 1.0f)) * s; // use s to express phi as -PI~PI
    float C = (phi - radians(restlength)) * s; // restlength represent degree here


    float dsum = 0.0, gamma = 1.0;
    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    // dont forget bend damping here
    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    float dL = (-C - alpha * l - dsum) / (gamma * gradC_invM_GradC + alpha);

    p0 += dL * invmass0 * gradC0;
    p1 += dL * invmass1 * gradC1;
    p2 += dL * invmass2 * gradC2;
    p3 += dL * invmass3 * gradC3;
    vstore3(p0, pt0, P);
    vstore3(p1, pt1, P);
    vstore3(p2, pt2, P);
    vstore3(p3, pt3, P);
    L[idx * 3] += dL;

}





// totally not understand what the SINGLE_WORL doing
kernel void constraintUpdate(
    int color_offset,
    int color_length, 

    float timeinc, 

    int type_length, 
    global int * type,
    
    int pts_length, 
    global int * pts_index, 
    global int * pts,
    
    int restlength_length, 
    global float * restlength ,
    
    int stiffness_length, 
    global float * stiffness ,
    
    int dampingratio_length, 
    global float * dampingratio ,
    
    int L_length, 
    global float * L ,
    
    int P_length, 
    global float * P ,
    
    int pprev_length, 
    global float * pprev ,
    
    int mass_length, 
    global float * mass 

) {

    if (timeinc == 0.0) return;

    int idx = get_global_id(0);
    if (idx >= color_length) return;
    idx += color_offset;
    
    int ctype = type[idx];
    int ptidx = pts_index[idx];
    float restlen = restlength[idx];
    float kstiff = stiffness[idx];
    float kdampratio = dampingratio[idx];

    if (kstiff == 0.0f) return;

    if (ctype == DISTANCE) {
        distanceUpdateXPBD(timeinc, idx, ptidx, pts, L, P, pprev, mass, restlen, kstiff, kdampratio);
    }

    if (ctype == BEND) {
        bendUpdateXPBD(timeinc, idx, ptidx, pts, L, P, pprev, mass, restlen, kstiff, kdampratio);
    }

}




