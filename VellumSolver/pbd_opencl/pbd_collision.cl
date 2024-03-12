
// #include </opt/hfs19.5/houdini/ocl/include/typedefines.h>
// #include </opt/hfs19.5/houdini/ocl/include/matrix.h>
// #include </opt/hfs19.5/houdini/ocl/include/quaternion.h>
// #include </opt/hfs19.5/houdini/ocl/include/reduce.h>


// this is so akward? chatgpt afterwards.
kernel void surfaceCollision( 
    float timeinc, 
    int P_length, 
    global float * P ,
    int pprevious_length, 
    global float * pprevious ,
    int mass_length, 
    global float * mass ,
    int hitnum_length, 
    global int * hitnum ,
    int hitpos_length, 
    global float * hitpos ,
    int hitnml_length, 
    global float * hitnml ,
    int hitv_length, 
    global float * hitv ,
    float  frictionscale 
) {
    int idx = get_global_id(0);
    if (idx >= P_length) return;

    if (hitnum[idx] <= 0 || mass[idx] == 0.0) return;

    float3 pos = vload3(idx, P);
    float3 pos_pre = vload3(idx, pprevious);
    float3 dp_hit = vload3(idx, hitv) * timeinc;
    float3 hit_normal = vload3(idx, hitnml);  

    // compute friction  
    float3 dp_prev = pos - pos_pre;
    float3 dp_prev_hit = dp_prev - dp_hit;
    float3 dp_normal = dot(dp_prev_hit, hit_normal) * hit_normal;
    float3 dp_tan = dp_prev_hit - dp_normal;
    float len_dp_tan = length(dp_tan);
    float len_dp_normal = length(dp_normal);

    // len_dp_normal is our approximate normal 
    // force, so we scale our tangent
    // by the ratio, clamping at 1.
    float fkin = 0.0;
    fkin = select((float)0, (float)1, (int)(fkin >= len_dp_tan));


}



kernel void projectCollision( 
    float timeinc, 
    int P_length, 
    global float * P ,
    int hitpos_length, 
    global float * hitpos ,
    int hitnml_length, 
    global float * hitnml ,
    int hitnum_length, 
    global int * hitnum ,
    int mass_length, 
    global float * mass 
) {
    int idx = get_global_id(0);
    if (idx >= P_length) return;

    if (hitnum[idx] <= 0 || mass[idx] == 0.0) return;

    float3 pos = vload3(idx, P);
    float3 hit_pos = vload3(idx, hitpos);
    float3 hit_normal = vload3(idx, hitnml);

    float C = dot(pos - hit_pos, hit_normal); // project to hitnormal and get intersection length
    pos -= min(C, 0.0f) * hit_normal; // throw the intersection outside
    vstore3(pos, idx, P);

}



