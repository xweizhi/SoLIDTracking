//SoLIDTracking
#include "SoLIDFieldMap.h"

SoLIDFieldMap * SoLIDFieldMap::fInstance = NULL;

//__________________________________________________________________
SoLIDFieldMap::SoLIDFieldMap()
{
    for (int i=0; i<ZSIZE; i++){
        for (int j=0; j<RSIZE; j++) { 
          Bz[i][j] = 0;
          Br[i][j] = 0;
        }
    }
  
    for (int i=0; i<ZSIZE_TARGET; i++){
        for (int j=0; j<RSIZE_TARGET; j++) { 
          Bz_target[i][j] = 0;
          Br_target[i][j] = 0;
        }
    }
    
    fField.SetXYZ(0.,0.,0.);
    targetFieldFlag = false;
  
    LoadFieldMap();
}
//__________________________________________________________________
void SoLIDFieldMap::LoadFieldMap()
{
  cout<<"loading SoLID field map"<<endl;
  ifstream infile;
  infile.open("solenoid_CLEOv8.dat");
  if (!infile.is_open()){
    cout<<"cannot open field map file"<<endl;
    exit(0);
  }
  double input[4];
  while(1){
    infile>>input[0]>>input[1]>>input[2]>>input[3];
    if (infile.eof()) break;
    
    int z_pos = fabs( (int)(input[1] + ZSHIFT) );
    int r_pos = fabs( (int)(input[0]) );
    Bz[z_pos][r_pos] = input[3]/1000.;//convert from gauss to kgauss
    Br[z_pos][r_pos] = input[2]/1000.;//convert from gauss to kgauss
  }
  
  infile.close();
  
}
//__________________________________________________________________
void SoLIDFieldMap::LoadTargetFieldMap()
{
    cout<<"loading SoLID target field map"<<endl;
    ifstream infile;
    infile.open("oxford_ptarget.dat");
    if (!infile.is_open()){
        cout<<"cannot open target field map file"<<endl;
        exit(0);
    }
    double input[4];
    while(1){
        infile>>input[0]>>input[1]>>input[2]>>input[3];
        if (infile.eof()) break;

        int z_pos = fabs( (int)(input[1] + ZSHIFT_TARGET) );
        int r_pos = fabs( (int)(input[0]) );
        Bz_target[z_pos][r_pos] = input[3]/1000.;//convert from gauss to kgauss
        Br_target[z_pos][r_pos] = input[2]/1000.;//convert from gauss to kgauss
    }

    infile.close();
    targetFieldFlag = true;
  
}
//___________________________________________________________________
TVector3 & SoLIDFieldMap::GetBField(double x, double y, double z)
{
    //here use cm, other place use m
    x = 100*x;
    y = 100*y;
    z = 100*z + ZSHIFT;
    double r = sqrt(x*x + y*y);
    if (r >= RSIZE || z <= 0 || z >= ZSIZE){
        fField.Clear();
        return fField;
    }else{

        int z_max, z_min, r_max, r_min;
        r_min = (int)r;
        r_max = r_min + RSTEP;
        z_min = (int)z;
        z_max = z_min + ZSTEP;

        double f21_Bz =  Bz[z_max][r_min];
        double f21_Br =  Br[z_max][r_min];
        double f22_Bz =  Bz[z_max][r_max];
        double f22_Br =  Br[z_max][r_max];

        double f11_Bz =  Bz[z_min][r_min];
        double f11_Br =  Br[z_min][r_min];
        double f12_Bz =  Bz[z_min][r_max];
        double f12_Br =  Br[z_min][r_max];

        //linear interpolation
        double Bzi = (1./((r_max - r_min)*(z_max - z_min)))*(f11_Bz*( z_max - z )*( r_max - r ) +
                                                             f21_Bz*( z - z_min )*( r_max - r ) +
                                                             f12_Bz*( z_max - z )*( r - r_min ) +
                                                             f22_Bz*( z - z_min )*( r - r_min ) );

        double Bri = (1./((r_max - r_min)*(z_max - z_min)))*(f11_Br*( z_max - z )*( r_max - r ) +
                                                             f21_Br*( z - z_min )*( r_max - r ) +
                                                             f12_Br*( z_max - z )*( r - r_min ) +
                                                             f22_Br*( z - z_min )*( r - r_min ) );

        fField.SetXYZ(Bri*x/r, Bri*y/r, Bzi);

        if (targetFieldFlag){
            //super impose the target field map
            z = z - (ZSHIFT - 350);
            r = sqrt(z*z + y*y);
            x = x + 200;
            
            if ( r < 200 && x >=0 && x < 401 && fabs(z) < 65 && fabs(y) < 65 && fabs(x - 200) < 65){
                fField.Clear();
                int x_min = (int)x;
                int x_max = x_min + ZSTEP_TARGET;
                r_min = (int)r;
                r_max = r_min + RSTEP_TARGET;
                f21_Bz =  Bz_target[x_max][r_min];
                f21_Br =  Br_target[x_max][r_min];
                f22_Bz =  Bz_target[x_max][r_max];
                f22_Br =  Br_target[x_max][r_max];

                f11_Bz =  Bz_target[x_min][r_min];
                f11_Br =  Br_target[x_min][r_min];
                f12_Bz =  Bz_target[x_min][r_max];
                f12_Br =  Br_target[x_min][r_max];
                
                Bzi = (1./((r_max - r_min)*(x_max - x_min)))*(f11_Bz*( x_max - x )*( r_max - r ) +
                                                                 f21_Bz*( x - x_min )*( r_max - r ) +
                                                                 f12_Bz*( x_max - x )*( r - r_min ) +
                                                                 f22_Bz*( x - x_min )*( r - r_min ) );

                Bri = (1./((r_max - r_min)*(x_max - x_min)))*(f11_Br*( x_max - x )*( r_max - r ) +
                                                                 f21_Br*( x - x_min )*( r_max - r ) +
                                                                 f12_Br*( x_max - x )*( r - r_min ) +
                                                                 f22_Br*( x - x_min )*( r - r_min ) );
                
                //so the z-axis of target field is pointing along the x-axis in solid coordinate
                //fField.SetXYZ(fField.X() + Bri*x/r, fField.Y() + Bri*y/r, fField.Z() + Bzi); 
                fField.SetXYZ(fField.X() + Bzi, fField.Y() + Bri*y/r, fField.Z() + Bri*z/r); 
            }
        }

        return fField;
    }

}














