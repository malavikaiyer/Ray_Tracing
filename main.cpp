#include<iostream>
#include<stdio.h>
#include<math.h>
#include<GLUT/glut.h>


using namespace std;

#define MAX_SIZE 512

struct	points
{
    float   x, y, z;
};

typedef struct	rgb_struct
{
    float   r, g, b;
} rgb;

struct Plane
{
    float a,b,c,d;
    float x_min, y_min, z_min, x_max, y_max, z_max;
    float refr_index;
};

struct Sphere{
    float cx,cy,cz;
    float r;
    float refr_index;
};

struct	points	From, At, up;
float	VXR, VXL, VYB, VYT;
float	ax, ay, az, bx, by, bz, cx, cy, cz;
float	viewangle, angle, tanv2;
float	xinterval, yinterval;


float lsx, lsy, lsz;
rgb	il, ia;
rgb	ka1, kd1, ks1;
rgb	ka2, kd2, ks2;
rgb	ka3, kd3, ks3;
rgb	ka4, kd4, ks4;
rgb	ka5, kd5, ks5;
rgb	ka6, kd6, ks6;
rgb	ka7, kd7, ks7;
rgb	tka1, tkd1, tka2, tkd2, tka3, tkd3;
rgb	tka4, tkd4, tka5, tkd5, tka6, tkd6;
rgb tka7, tkd7;

int	phong1, phong2, phong3, phong4, phong5, phong6;


int		xmax_pixel, ymax_pixel;


float *texture_R;
float *texture_G;
float *texture_B;
float noise_tabl[65][65];

Plane p_arr[20];
int numPlane = 0;
Sphere s_arr[10];
int numSphere = 0;

int TABLE_SIZE_1 = 60;
int TABLE_SIZE = 60;

FILE	*outpfile;


int Read_Information()
{
    string str;
    
    if ( freopen("input.dat","r", stdin) == NULL)
    {
        printf("ERROR: Could not open input.dat for read!\n");
        return(0);
    }
    
    while(cin>>str){
        
        if(str == "\n")continue;
        
        else if(str == "From" || str == "from")
        {
            cin>>From.x>>From.y>>From.z;
        }
        else if(str == "At" || str == "at")
        {
            cin>>At.x>>At.y>>At.z;
        }
        else if(str == "Up" || str == "up")
        {
            cin>>up.x>>up.y>>up.z;
        }
        else if(str == "ViewAngle" || str == "viewangle")
        {
            cin>>viewangle;
            angle = viewangle * 3.14/180.0;
            tanv2 = tan(angle/2.0);
        }
        else if(str == "Viewport" || str == "viewport")
        {
            cin>>VXL>>VXR>>VYB>>VYT;
        }
        else if(str == "Light" || str == "light")
        {
            cin>>lsx>>lsy>>lsz;
        }
        else if(str == "Plane" || str == "plane")
        {
            cin>>p_arr[numPlane].a>>p_arr[numPlane].b>>p_arr[numPlane].c;
            cin>>p_arr[numPlane].d;
            cin>>p_arr[numPlane].x_min>>p_arr[numPlane].y_min>>p_arr[numPlane].z_min;
            cin>>p_arr[numPlane].x_max>>p_arr[numPlane].y_max>>p_arr[numPlane].z_max;
            cin>>p_arr[numPlane].refr_index;
            numPlane++;
        }
        else if(str == "Sphere" || str == "sphere")
        {
            cin>>s_arr[numSphere].cx>>s_arr[numSphere].cy>>s_arr[numSphere].cz>>s_arr[numSphere].r;
            cin>>s_arr[numSphere].refr_index;
            numSphere++;
        }
        else if(str == "ImageSize" || str == "imagesize")
        {
            cin>>xmax_pixel>>ymax_pixel;
            if (xmax_pixel > MAX_SIZE || ymax_pixel > MAX_SIZE)
            {
                printf("Error: Exceeded max image size %d x %d\n", xmax_pixel, ymax_pixel);
                printf("Reset to max image size: %d x %d\n", MAX_SIZE, MAX_SIZE);
                xmax_pixel = MAX_SIZE-1;
                ymax_pixel = MAX_SIZE - 1;
            }
        }
    }
    fclose(stdin);
    
    if ((outpfile = fopen("output.out","wb")) == NULL) {
        printf("ERROR:  cannot open output.out for write.\n");
        return(0);
    }
    
    texture_R = new float [xmax_pixel * ymax_pixel ];
    texture_G = new float [xmax_pixel * ymax_pixel ];
    texture_B = new float [xmax_pixel * ymax_pixel ];
    
    
    return(1);
}



void Normalize(float *x,float *y,float *z)
{
    float	norm;
    
    norm = sqrt( *x * *x + *y * *y + *z * *z );
    if (norm != 0.0) {
        *x = *x / norm;
        *y = *y / norm;
        *z = *z / norm;
    }
}


float 	Power(float base,int exp)
{
    int	i;
    float	value;
    
    value = 1.0;
    for (i=1; i<=exp; i++)
        value *= base;
    
    return( value );
}


void Compute_M()
{
    
    cx = At.x - From.x;
    cy = At.y - From.y;
    cz = At.z - From.z;
    Normalize(&cx, &cy, &cz);
    
    
    ax = cy*up.z - up.y*cz;
    ay = up.x*cz - cx*up.z;
    az = cx*up.y - up.x*cy;
    Normalize(&ax, &ay, &az);
    
    
    bx = ay*cz - cy*az;
    by = cx*az - ax*cz;
    bz = ax*cy - cx*ay;
}


void Setup_Parameters()
{
    
    
    
    Compute_M();
    
    Normalize(&lsx, &lsy, &lsz);
    
    xinterval = (VXR - VXL) / xmax_pixel;
    yinterval = (VYT - VYB) / ymax_pixel;
    
    
    il.r = 1.0;	il.g = 1.0;	il.b = 1.0;
    ia.r = 1.0;	ia.g = 1.0;	ia.b = 1.0;
    
    
    ka1.r = 0.3;	ka1.g = 0.0;	ka1.b = 0.0;
    kd1.r = 0.7;	kd1.g = 0.0;	kd1.b = 0.0;
    ks1.r = 1.0;	ks1.g = 1.0;	ks1.b = 1.0;
    tka1.r = 0.2;	tka1.g = 0.0;	tka1.b = 0.0;
    tkd1.r = 0.2;	tkd1.g = 0.0;	tkd1.b = 0.0;
    phong1 = 60;
    
    ka2.r = 0.0;	ka2.g = 0.3;	ka2.b = 0.0;
    kd2.r = 0.0;	kd2.g = 0.7;	kd2.b = 0.0;
    ks2.r = 1.0;	ks2.g = 1.0;	ks2.b = 1.0;
    tka2.r = 0.0;	tka2.g = 0.2;	tka2.b = 0.0;
    tkd2.r = 0.0;	tkd2.g = 0.2;	tkd2.b = 0.0;
    phong2 = 90;
    
    ka3.r = 0.0;	ka3.g = 0.0;	ka3.b = 0.3;
    kd3.r = 0.0;	kd3.g = 0.0;	kd3.b = 0.7;
    ks3.r = 1.0;	ks3.g = 1.0;	ks3.b = 1.0;
    tka3.r = 0.0;	tka3.g = 0.0;	tka3.b = 0.2;
    tkd3.r = 0.0;	tkd3.g = 0.0;	tkd3.b = 0.2;
    phong3 = 120;
    
    
    ka4.r = 0.1;	ka4.g = 0.1;	ka4.b = 0.1;
    kd4.r = 0.2;	kd4.g = 0.5;	kd4.b = 0.0;
    ks4.r = 1.0;	ks4.g = 1.0;	ks4.b = 1.0;
    tka4.r = 0.0;	tka4.g = 0.1;	tka4.b = 0.0;
    tkd4.r = 0.0;	tkd4.g = 0.7;	tkd4.b = 0.0;
    phong4 = 120;
    
    ka5.r = 0.1;	ka5.g = 0.0;	ka5.b = 0.1;
    kd5.r = 0.7;	kd5.g = 0.0;	kd5.b = 0.7;
    ks5.r = 1.0;	ks5.g = 1.0;	ks5.b = 1.0;
    tka5.r = 0.0;	tka5.g = 0.0;	tka5.b = 0.1;
    tkd5.r = 0.0;	tkd5.g = 0.0;	tkd5.b = 0.7;
    phong5 = 120;
    
    
    ka6.r = 1.0;	ka6.g = 1.0;	ka6.b = 1.0;
    kd6.r = 0.8;	kd6.g = 0.8;	kd6.b = 0.8;
    ks6.r = 0.8;	ks6.g = 0.8;	ks6.b = 0.8;
    tka6.r = 0.0;	tka6.g = 0.0;	tka6.b = 0.5;
    tkd6.r = 0.0;	tkd6.g = 0.0;	tkd6.b = 0.0;
    phong6 = 120;
    
    
    ka7.r = 0.2;	ka7.g = 0.2;     ka7.b = 0.07;
    kd7.r = 0.7;	kd7.g = 0.6;	kd7.b = 0.22648;
    ks7.r = 0.62;	ks7.g = 0.55;	ks7.b = 0.36;
    tka7.r = 0.0;       tka7.g = 0.0;       tka7.b = 0.0;
    tkd7.r = 0.0;       tkd7.g = 0.0;       tkd7.b = 0.0;
}




float calc_noise(float iu, float iv, int direction)
{
    int i, j, x, y, left, right;
    float noise;
    float val;
    
    
    for (i = 0; i < 65; i++)
        for (j = 0; j < 65; j++)
        {
            val = rand() % 255;
            noise_tabl[i][j] = val / 256;
        }
    
    i = (int)iu;
    j = (int)iv;
    x = i%TABLE_SIZE;
    y = j%TABLE_SIZE;
    
    if (direction == 1)
    {
        if (x <= 0)
            left = 0;
        else
            left = x - 1;
        
        if (x >= TABLE_SIZE_1)
            right = TABLE_SIZE_1;
        else
            right = x + 1;
        
        noise = (noise_tabl[right][y] - noise_tabl[left][y]) / 2.0;
    }
    
    else {
        if (y <= 0)
            left = 0;
        else
            left = y - 1;
        
        if (y >= TABLE_SIZE_1)
            right = TABLE_SIZE_1;
        else
            right = y + 1;
        noise = (noise_tabl[x][right] - noise_tabl[x][left]) / 2.0;
    }
    return(noise);
    
}


void Bump_Map(float x, float y, float z, float xc, float yc, float zc, float r, float *nx, float *ny, float *nz, int methd)
{
    
    float xp, yp, zp, iu, iv, xu, yu, zu, xv, yv, zv;
    float fu, fv, a, dx, dy, dz, u, v, nnx, nny, nnz;
    
    xp = (x - xc) / r;
    yp = (y - yc) / r;
    zp = (z - zc) / r;
    
    u = asin(zp);
    v = atan2(yp, xp);
    
    iu = (u + 3.14) / (2 * 3.14) * (TABLE_SIZE_1);
    iv = (v + (3.14 / 2)) / 3.14 * (TABLE_SIZE_1);
    
    xu = -r*cos(v)*sin(u);
    xv = -r*sin(v)*cos(u);
    yu = r*cos(v)*cos(u);
    yv = -r*sin(v)*sin(u);
    zu = 0.0;
    zv = r*cos(v);
    
    fu = calc_noise(iu, iv, 1);
    fv = calc_noise(iu, iv, 2);
    
    if(methd == 1){
        dx = fu*xu + fv*xv;
        dy = fu*yu + fv*yv;
        dz = fu*zu + fv*zv;
    }
    
    else{
        nnx = *nx;
        nny = *ny;
        nnz = *nz;
        Normalize(&nnx, &nny, &nnz);
        
        
        dx = fv* (yu*nnz - nny*zu);
        dy = fv* (nnx*zu - xu*nnz);
        dz = fv* (xu*nny - nnx*yu);
        
        dx += fu* (yv*nnz - nny*zv);
        dy += fu* (nnx*zv - xv*nnz);
        dz += fu* (xv*nny - nnx*yv);
    }
    Normalize(&dx, &dy, &dz);
    
    a = sqrt(fu*fu + fv*fv);
    dx *= a;
    dy *= a;
    dz *= a;
    
    *nx += dx;
    *ny += dy;
    *nz += dz;
}



void Check_Sphere(float px,float py,float pz,float dx,float dy,float dz,float xc,float yc,float zc,float r,
                  float *t1, float *t2)
{
    float	a, b, c, xdiff, ydiff, zdiff, discr;
    
    xdiff = px-xc;
    ydiff = py-yc;
    zdiff = pz-zc;
    a = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff - r*r;
    b = 2.0*( dx*xdiff + dy*ydiff + dz*zdiff );
    c = dx*dx + dy*dy + dz*dz;
    
    
    discr = b*b - 4.0*a*c;
    if (discr < 0.0)
    {
        *t1 = -1.0;
        *t2 = -1.0;
    }
    else if (discr == 0.0)
    {
        *t1 = -b / (2.0*c);
        *t2 = -1.0;
    }
    else {
        discr = sqrt(discr);
        *t1 = (-b + discr) / (2.0*c);
        *t2 = (-b - discr) / (2.0*c);
    }
}



void Check_Plane(float px,float py,float pz,float dx,float dy,float dz,float a,float b,float c,float d,float *t1)
{
    *t1 = (-a*px - b*py - c*pz - d) / (a*dx + b*dy + c*dz);
}


int Check_Shadow(float px, float py, float pz)
{
    
    float t1, t2;
    
    for(int i = 0; i<numSphere; i++)
    {
        Check_Sphere(px, py, pz, lsx, lsy, lsz, s_arr[i].cx, s_arr[i].cy, s_arr[i].cz , s_arr[i].r, &t1, &t2);
        
        if(t1>= 0.0001 || t2>=0.0001)return 1;
    }
    return 0;
}


void Compute_Intersection(float px,float py,float pz,float dx,float dy, float dz,float t,float *newx,float *newy,float *newz)
{
    *newx = px + t*dx;
    *newy = py + t*dy;
    *newz = pz + t*dz;
}



void Compute_Color(int shadow_flag, float ipx,float ipy,float  ipz,float  nx,float  ny,float  nz,
                   rgb ia,rgb ka,rgb kd, rgb ks,int n,float *r,float *g, float *b)
{
    float	vx, vy, vz, rx, ry, rz;
    float	ndotl, vdotr, cosalphapower;
    
    
    vx = From.x - ipx;
    vy = From.y - ipy;
    vz = From.z - ipz;
    Normalize(&vx, &vy, &vz);
    
    ndotl = nx*lsx + ny*lsy + nz*lsz;
    rx = 2.0*ndotl*nx - lsx;
    ry = 2.0*ndotl*ny - lsy;
    rz = 2.0*ndotl*nz - lsz;
    
    vdotr = vx*rx + vy*ry + vz*rz;
    
    *r = ia.r * ka.r;
    *g = ia.g * ka.g;
    *b = ia.b * ka.b;
    
    
    if (ndotl >= 0.0 && shadow_flag==0)
    {
        *r = *r + kd.r*ndotl*il.r;
        *g = *g + kd.g*ndotl*il.g;
        *b = *b + kd.b*ndotl*il.b;
        
        if (vdotr >= 0.0)
        {
            cosalphapower = Power(vdotr, n);
            
            *r = *r + ks.r*cosalphapower*il.r;
            *g = *g + ks.g*cosalphapower*il.g;
            *b = *b + ks.b*cosalphapower*il.b;
        }
    }
    
    if (*r > 1.0) *r = 1.0;
    if (*g > 1.0) *g = 1.0;
    if (*b > 1.0) *b = 1.0;
}


void Ray_Reflection(float px, float py, float pz, float dx, float dy, float dz, float nx, float ny, float nz, float *rlx, float *rly, float *rlz)
{
    float ix = dx - px;
    float iy = dy - py;
    float iz = dz - pz;
    float c1 = -(nx*dx + ny*dy + nz*dz);
    *rlx = dx + ( 2 * nx * c1);
    *rly = dy + ( 2 * ny * c1);
    *rlz = dz + ( 2 * nz * c1);
    
}


void Ray_Refraction(float dx, float dy, float dz, float nx, float ny, float nz, float *rrx, float *rry, float *rrz, float n2, float n1)
{
    
    //n1-> original medium
    //n2-> new medium
    float n, c1, c2;
    n = n1 / n2;
    c1 = - ((nx*dx)+(ny*dy)+(nz*dz));
    c2 = sqrt( 1 - (n * n) * (1 - (c1 * c1)));
    *rrx = (n * dx) + (n * c1 - c2) * nx;
    *rry = (n * dy) + (n * c1 - c2) * ny;
    *rrz = (n * dz) + (n * c1 - c2) * nz;
}


rgb Recursive_Ray_Trace(int level, points frm, float dx, float dy, float dz, float refr_index, int obj_type, int obj_id)
{
    
    int	    obj, obj_num, shadow_flag;
    int	    texture;
    float	nx, ny, nz;
    float	t_min, t1, t2, ipx, ipy, ipz;
    float   rlx, rly, rlz, rrx, rry, rrz;
    float	r, g, b;
    int     i;
    rgb     newcolor, reflect_color, refract_color;
    float   refraction_index = refr_index;
    
    t_min = 999.0;
    obj_num = 0;
    obj = 0;
    texture = 0;
    
    //For intersection with spheres
    for(i = 0; i<numSphere; i++)
    {
        if(obj_type == 1 && obj_id == i)continue;
        Check_Sphere(frm.x, frm.y, frm.z, dx, dy, dz, s_arr[i].cx, s_arr[i].cy, s_arr[i].cz ,s_arr[i].r, &t1, &t2);
        
        if (t1>=0.0)
        {
            t_min = t1;
            obj = 1;
            obj_num = i;
            refraction_index = s_arr[i].refr_index;
            Compute_Intersection(frm.x, frm.y, frm.z, dx, dy, dz, t1, &ipx, &ipy, &ipz);
        }
        
        if (t2>=0.0 && t2<t_min)
        {
            t_min = t2;
            obj = 1;
            obj_num = i;
            refraction_index = s_arr[i].refr_index;
            Compute_Intersection(frm.x, frm.y, frm.z, dx, dy, dz, t2, &ipx, &ipy, &ipz);
        }
        
    }
    
    //For intersection with planes
    for(i=0; i<numPlane; i++)
    {
        if(obj_type == 2 && obj_id == i)
            continue;
        Check_Plane(frm.x, frm.y, frm.z, dx, dy, dz, p_arr[i].a, p_arr[i].b, p_arr[i].c, p_arr[i].d, &t1);
        
        if (t1 >= 0.0 && t1<t_min)
        {
            Compute_Intersection(frm.x, frm.y, frm.z, dx, dy, dz, t1, &ipx, &ipy, &ipz);
            
            if (ipx >= p_arr[i].x_min && ipx <= p_arr[i].x_max &&
                ipy >= p_arr[i].y_min  && ipy <= p_arr[i].y_max &&
                ipz >=  p_arr[i].z_min && ipz <= p_arr[i].z_max )
            {
                
                t_min = t1;
                obj_num = i;
                obj = 2;
                refraction_index = p_arr[i].refr_index;
            }
        }
    }
    
    
    if(fabs(frm.x - ipx) < 0.00001 && fabs(frm.y - ipy) < 0.00001 && fabs(frm.z - ipz) < 0.00001)
        obj = 0;
    
    switch (obj)
    {
        case 0 :
            r = 0.3;
            g = 0.3;
            b = 0.5;
            break;
            
        case 1 :
            //Sphere
            nx = ipx - s_arr[obj_num].cx;
            ny = ipy - s_arr[obj_num].cy;
            nz = ipz - s_arr[obj_num].cz;
            Normalize(&nx, &ny, &nz);
            
            shadow_flag = 0;
            shadow_flag = Check_Shadow(ipx, ipy, ipz );
            texture = 0;
            
            if(obj_num == 0)
            {
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka1, kd1, ks1, phong1, &r, &g, &b);
            }
            if(obj_num == 1)
            {
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka2, kd2, ks2, phong2, &r, &g, &b);
            }
            if(obj_num == 2)
            {
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka3, kd3, ks3, phong3, &r, &g, &b);
            }
            if(obj_num == 3)
            {
                if(rand()%3 == 1)
                    Bump_Map(ipx, ipy, ipz, dx, dy, dz, 1, &nx, &ny, &nz,1);
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka7, kd7, ks7, phong6, &r, &g, &b);
            }
            break;
            
            
        case 2 :
            //Planes
            nx = p_arr[obj_num].a;
            ny = p_arr[obj_num].b;
            nz = p_arr[obj_num].c;
            
            shadow_flag = 0;
            shadow_flag = Check_Shadow( ipx, ipy, ipz );
            
            if(obj_num == 0)
            {
                
                if((((int)ipx % 4) < 2 && ((int)ipz % 4) < 2)|| (((int)ipx % 4) >= 2 && ((int)ipz % 4) >= 2)){
                    texture = 0;
                }
                else texture = 1;
                
                if (texture == 1)
                {
                    Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka5, tkd5, ks5, phong6, &r, &g, &b);
                }
                else{
                    Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka6, kd6, ks6, phong6, &r, &g, &b);
                }
            }
            if(obj_num == 1)
            {
                //Wood grain
                float radius1, radius2, ang;
                float u, v, w;
                int grain;
                u = (p_arr[1].x_max - p_arr[1].x_min)/2 + ipx;
                v = (p_arr[1].y_max - p_arr[1].y_min)/2 + ipy;
                w = (p_arr[1].z_max - p_arr[1].z_min)/2 + ipz;
                radius1 = sqrt(u*u + v * v);
                ang = atan(u/v);
                radius2 = radius1 + 2 * sin(20 * ang) + u / 15;
                grain = (int)radius2 % 4;
                if(grain <= 3 )
                    Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka4, kd4, ks4, phong4, &r, &g, &b);
                else
                    Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka6, tkd6, ks6, phong4, &r, &g, &b);
                
            }
            if(obj_num  == 2)
            {
                if(rand() % 2 == 1)
                {
                    Bump_Map(ipx, ipy, ipz, dx, dy, dz, 1, &nx, &ny, &nz,1);
                    Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka3, kd3, ks3, phong3, &r, &g, &b);
                }
                else
                    Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka3, kd3, ks3, phong3, &r, &g, &b);
            }
            
    }
    
    
    Ray_Reflection(ipx,ipy,ipz, dx, dy, dz, nx, ny, nz, &rlx, &rly, &rlz);
    if(refr_index > 1)
        Ray_Refraction(dx, dy, dz, nx, ny, nz, &rrx, &rry, &rrz, refr_index, 1);
    else
        Ray_Refraction(dx, dy, dz, nx, ny, nz, &rrx, &rry, &rrz, refr_index, refraction_index);
    
    points newFrom;
    newFrom.x = ipx;
    newFrom.y = ipy;
    newFrom.z = ipz;
    
    if(level >= 4 || obj == 0)
    {
        newcolor.r = r ;
        newcolor.g = g ;
        newcolor.b = b ;
        return newcolor;
    }
    level += 1;
    
    
    if( obj == 1 && obj_num == 2)
    {
        refract_color = Recursive_Ray_Trace(level, newFrom, rrx, rry, rrz, refraction_index,obj, obj_num);
        newcolor.r = r*0.1 +  refract_color.r/(1.1*level);
        newcolor.g = g*0.1 +  refract_color.g/(1.1*level);
        newcolor.b = b*0.1 +  refract_color.r/(1.1*level);
    }
    else if(obj == 2 && obj_num == 2)
    {
        newcolor.r = r ;
        newcolor.g = g ;
        newcolor.b = b;
    }
    else
    {
        reflect_color = Recursive_Ray_Trace(level, newFrom, rlx, rly, rlz, 1, obj, obj_num);
        newcolor.r = r + reflect_color.r/(1.2*level) ;
        newcolor.g = g + reflect_color.g/(1.2*level);
        newcolor.b = b + reflect_color.b/(1.2*level);
    }
    return newcolor;
}


void Ray_Generate(){
    
    int	    xp, yp;
    int	    buf_ptr;
    float	xv, yv, dx, dy, dz;
    float   u, v;
    rgb     newcolor;
    
    buf_ptr = 0;
    for (xp=0; xp<xmax_pixel; xp++)
    {
        u = (float)xp/xmax_pixel;
        
        for (yp=0; yp<ymax_pixel; yp++)
        {
            v = (float)yp/ymax_pixel;
            
            xv = VXL + xp * xinterval;
            yv = VYB + yp * yinterval;
            
            
            dx = ax*xv*tanv2 + bx*yv*tanv2 + cx;
            dy = ay*xv*tanv2 + by*yv*tanv2 + cy;
            dz = az*xv*tanv2 + bz*yv*tanv2 + cz;
            
            
            newcolor = Recursive_Ray_Trace(1, From, dx, dy, dz, 1, -1, -1);
            
            texture_R[xp + xmax_pixel * yp] = newcolor.r;
            texture_G[xp + xmax_pixel * yp] = newcolor.g;
            texture_B[xp + xmax_pixel * yp] = newcolor.b;
            
        }
    }
    
    printf("Writing to image...\n");
    fwrite(&xmax_pixel, sizeof(int), 1, outpfile);
    fwrite(&ymax_pixel, sizeof(int), 1, outpfile);
    
    fwrite(texture_R, sizeof(float), xmax_pixel*ymax_pixel, outpfile);
    fwrite(texture_G, sizeof(float), xmax_pixel*ymax_pixel, outpfile);
    fwrite(texture_B, sizeof(float), xmax_pixel*ymax_pixel, outpfile);
    
    fclose(outpfile);
}


void myinit(void)
{
    
    
    glClearColor(0.0, 0.0, 0.0, 1.0);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 512.0, 0.0, 512.0);
    glMatrixMode(GL_MODELVIEW);
}

void display( void )
{
    int s, t;
    float  r, g, b;
    
    glClear(GL_COLOR_BUFFER_BIT);
    
    for(t = 0; t < ymax_pixel; t++) {
        for(s = 0; s < xmax_pixel; s++) {
            
            r = texture_R[s + xmax_pixel * t];
            g = texture_G[s + xmax_pixel * t];
            b = texture_B[s + xmax_pixel * t];
            
            glColor3f(r, g, b);
            glBegin(GL_POINTS);
            glVertex2f(s,t);
            glEnd();
        }
    }
    
    glFlush();
}


int main(int argc, char**argv)
{
    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(500,500);
    glutInitWindowPosition(0,0);
    glutCreateWindow("Ray Tracing");
    glutDisplayFunc(display);
    
    Read_Information();
    Setup_Parameters();
    Ray_Generate();
    myinit();
    glutMainLoop();
    
    return(0);
}