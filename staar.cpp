///////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////   STAAR Software  //////////////////////////////////////////////
/////////////////  Initiated by Chandra Narasimhan  for Dr. Liz Howell ///////////////////
/////////////////    This software should not be distributed as open source //////////////
//////////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <locale>
#include <cstdio>
#include <cmath>
#include <unistd.h>
#include <getopt.h>

using namespace std;

#define LINE_LENGTH 100000
#define MAX_RESULTS 10000
#define MAX_FILE_NAME 1000
#define MAX_HET_ATOM 100000
#define MAX_ALT_LOC 10  

class Options{

public:
  bool   rCode;
  char currentDirectory[MAX_FILE_NAME];
  char *pdbfolder;
  char *pdbListfile;
  char *opFile;

  Options(){
    rCode = false;
    getcwd(currentDirectory, MAX_FILE_NAME);
    pdbfolder = NULL;
    pdbListfile = NULL;
    opFile = NULL;
  }

  Options(int argc, char **argv){
    rCode = false;
    getcwd( currentDirectory, MAX_FILE_NAME );
    pdbfolder = NULL;
    pdbListfile = NULL;
    opFile = NULL;
    if ( parseCmdline( argc, argv ) ){
      exit(1);
    }
  }
  
  void printHelp(){
    cerr << "-h or --help          " << "Displays this message" << endl;
    cerr << "-p or --pdblist       " << "Specifies the PDBList file" << endl;
    cerr << "-f or --pdbfolder     " << "Specifies the folder for PDB files" << endl;
    cerr << "-o or --out           " << "Specifies the output file" << endl;
    cerr << "-c or --center        " << "Specifies whether to calculat center of charge" << endl;
  }

  int parseCmdline(int argc, char **argv)
  {
  
    int c;
    static struct option long_options[] = 
      {
	{"help",      required_argument, 0, 'h'},
	{"pdblist",   required_argument, 0, 'p'},
	{"pdbfolder", required_argument, 0, 'f'},
	{"out",       required_argument, 0, 'o'}, 
	{"center",    required_argument, 0, 'c'},
	{0, 0, 0, 0}
      };
    int option_index;
    while( !( ( c = getopt_long(argc, argv, "hp:f:o:c", long_options, &option_index) ) < 0 ) )
      {
	switch(c)
	  {
	  case 'h':
	    // print some help info
	    printHelp();
	    exit(0);
	    break;

	  case 'p':
	    this->pdbListfile = optarg;
	    break;

	  case 'f':
	    this->pdbfolder = optarg;
	    break;

	  case 'o':
	    this->opFile = optarg;
	    break;

	  case 'c':
	    this->rCode = true;
	    break;
	  }
      }

    if( optind < argc ){
      cerr << "There are some extra commands inputted. Please cleanup the command line!" << endl;
      printHelp();
      // print some help info
      return -1;
    }else if( !pdbListfile ){
      cerr << "Must specify the PDB list file with -p or --pdblist" <<  endl;
      printHelp();
      return -1;
    }else if( !pdbfolder ){
      cerr << "Must specify the PDB folder with -f or --pdbfolder" <<  endl;
      printHelp();
      return -1;
    }else if( !opFile ){
      cerr << "Must specify the op file with -o or --op" <<  endl;
      printHelp();
      return -1;
    }

    return 0;

  }

};

//char rCode='N';  //// no center of charge calculation for carboxyl group
//char rCode='C';  //// center of charge calculation for carboxyl group
bool rCode;
////////-------------------------centroid calculation Constants-------------
float massC=12.0107,massN=14.00674,massO=15.9994;
float chargeC=0.898571,chargeO=-0.83462,chargeH=-0.229336;
int method=0;
////////----------------------String Parse Functions------------------------

char tmpSequence[100000];

int charToInt(char *ptr,int len){
  char line[200];
  memcpy(&line[0],ptr,len);
  line[len]='\0';
  return strtol(&line[0], (char **)NULL, 10);
}

float charToFloat(char *ptr,int len){
  char line[200];
  memcpy(&line[0],ptr,len);
  line[len]='\0';
  return (float) strtod(&line[0], (char **)NULL);
}

////////--------------------Analytical Geometry-----------------------------------

// Distance between 2 points in space

float findDistance(float x1,float y1,float z1,float x2,float y2,float z2){
  float sum;

  sum = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1)) ;
  return sum;
}

// The following function calculates the length of a vector

float getNormalValue(float x,float y,float z){
  return sqrt((x*x + y*y + z*z));
}

// The following function calculates the dot product between 2 vectors 

float getDotProduct(float x1,float y1,float z1,float x2,float y2,float z2){
  return (x1*x2 + y1*y2 + z1*z2);
}

class coord{
public:
  float x;
  float y;
  float z;
  
  coord(){
    x=0.0;
    y=0.0;
    z=0.0;    
  }
  
  coord(float x1,float y1,float z1){
    x=x1;
    y=y1;
    z=z1;
  }

  void setCoord(float x1,float y1,float z1){
    this->x=x1;
    this->y=y1;
    this->z=z1;
  }

  ~coord(){
  }
};

float getDeterminant(coord *a,coord *b,coord *c){
  float a1,a2,a3,b1,b2,b3,c1,c2,c3;
  float deter;

  a1=a->x; a2=a->y; a3=a->z;
  b1=b->x; b2=b->y; b3=b->z;
  c1=c->x; c2=c->y; c3=c->z;

  deter= a1*(b2*c3-b3*c2)-a2*(b1*c3-b3*c1)+a3*(b1*c2-b2*c1);

  return deter;
}

float getPlaneEquation(coord *p1,coord *p2,coord *p3,coord *result){
  float det,detX,detY,detZ;
  int error=1;
  coord a,b,c;

  a.setCoord(p1->x,p1->y,p1->z);
  b.setCoord(p2->x,p2->y,p2->z);
  c.setCoord(p3->x,p3->y,p3->z);

  det=getDeterminant(&a,&b,&c);
  a.setCoord(1.0,p1->y,p1->z);
  b.setCoord(1.0,p2->y,p2->z);
  c.setCoord(1.0,p3->y,p3->z);
  detX=getDeterminant(&a,&b,&c);

  a.setCoord(p1->x,1.0,p1->z);
  b.setCoord(p2->x,1.0,p2->z);
  c.setCoord(p3->x,1.0,p3->z);
  detY=getDeterminant(&a,&b,&c);

  a.setCoord(p1->x,p1->y,1.0);
  b.setCoord(p2->x,p2->y,1.0);
  c.setCoord(p3->x,p3->y,1.0);
  detZ=getDeterminant(&a,&b,&c);
    
  result->x=detX;
  result->y=detY;
  result->z=detZ;

  return det;
}

float getAngleBetweenPlaneAndLine(coord *plane,coord *pointX,coord *pointY){
  float normalV,normalP,dotProd,angle;
  float vec1_x,vec1_y,vec1_z;
  double cosValue;

  normalV=getNormalValue(plane->x,plane->y,plane->z);
  vec1_x=pointY->x-pointX->x;
  vec1_y=pointY->y-pointX->y;
  vec1_z=pointY->z-pointX->z;

  normalP=getNormalValue(vec1_x,vec1_y,vec1_z);
  dotProd=getDotProduct(plane->x,plane->y,plane->z,vec1_x,vec1_y,vec1_z);
  cosValue=dotProd/(normalV*normalP);
  if((cosValue<-1.00000)||(cosValue>1.0000000)) {
    cout <<" Error: Can't calculate cosine inverse. Cosine of angle lies outside -1 and 1 \n"<<flush;
    exit(20);
  }else{
    angle=acos(cosValue);
    angle=angle*180/3.1412;
    angle=90.0-angle;
    if(angle<0.0) angle=(-1.0)*angle;
  }
  return angle;
}

/// The following line calculates angle between 2 planes

float getAngleBetweenPlaneAndPlane(coord *plane1,coord *plane2){
  float normalV,normalP,dotProd,angle;
  double cosValue;

  normalV=getNormalValue(plane1->x,plane1->y,plane1->z);
  normalP=getNormalValue(plane2->x,plane2->y,plane2->z);

  dotProd=getDotProduct(plane1->x,plane1->y,plane1->z,plane2->x,plane2->y,plane2->z);
  cosValue=dotProd/(normalV*normalP);
  if((cosValue<-1.00000)||(cosValue>1.0000000)) {
    cout <<" Error: Can't calculate cosine inverse. Cosine of angle lies outside -1 and 1 \n"<<flush;
    exit(20);
  }else{
    angle=acos(cosValue);
    angle=angle*180/3.1412;
    if(angle>90.0) angle=180.0-angle;
  }
  return angle;
}

/// The following line calculates angle between 2 planes

float getConstantForPlaneLineIntercept(coord *plane,coord *point,float intercept){
  float t=0.0,a,b,c,x1,y1,z1,denom;

  a=plane->x;
  b=plane->y;
  c=plane->z;
  x1=point->x;
  y1=point->y;
  z1=point->z;

  denom=a*a+b*b+c*c;
  if(denom!=0.0){
    t=(-1.0)*(intercept+a*x1+b*y1+c*z1)/denom;
  }else{
    cout<< "Warning: Plane vector does not exist\n";
  }
  return t;
}

void getPlaneProjectCoordinate(coord *plane,coord *point,float intercept,coord *finalAns){
  float t,a,b,c,x1,y1,z1;

  t=getConstantForPlaneLineIntercept(plane,point,intercept);
  a=plane->x;
  b=plane->y;
  c=plane->z; 
  x1=point->x;
  y1=point->y;
  z1=point->z;
  finalAns->x=x1+a*t;
  finalAns->y=y1+b*t;
  finalAns->z=z1+c*t;

}

// The following function calculates the angle between 2 lines with intercept at (xc,yc,zc)
// (x1,y1,z1) should be a point in the plane of aromatic ring. (x2,y2,z2) is the coordinate of the
// interacting atom of interest. It returns angle in degrees.

float findAngle(float x1,float y1,float z1,float xc,float yc,float zc,float x2,float y2,float z2){
  float vec1_x,vec1_y,vec1_z,vec2_x,vec2_y,vec2_z;
  float lenVec1,lenVec2,prod,angle;
  double cosValue;
  
  vec1_x=x1-xc;
  vec1_y=y1-yc;
  vec1_z=z1-zc;  

  vec2_x=x2-xc;
  vec2_y=y2-yc;
  vec2_z=z2-zc;  

  lenVec1=getNormalValue(vec1_x,vec1_y,vec1_z);
  lenVec2=getNormalValue(vec2_x,vec2_y,vec2_z);
  prod=getDotProduct(vec1_x,vec1_y,vec1_z,vec2_x,vec2_y,vec2_z);

  cosValue=prod/(lenVec1*lenVec2);
  if((cosValue<-1.00000)||(cosValue>1.0000000)) {
    cout <<" Error: Can't calculate cosine inverse. Cosine of angle lies outside -1 and 1 \n"<<flush;
    angle=1000.0;
  }else{
    angle=acos(cosValue);
    angle=angle*180/3.1412;
    ///    if(angle>90.0) angle=180.0-angle;
  }
  
  return angle;
}





//------------------------------------------------------------------------------------------------------
// class amino acid stores information about all the atoms in the amino acid. The maximum total number of
// atoms in the amino acid is assumed to be 25. The field atomName stores the names of the atom. For each 
// atom, the number of charecters allocated is 4. Hence, the size of this array is 4*25=100.
// numberOfAtoms is a count of total atom information present in an aminoacid (this will be <25).
// 
// aminoAcid class stores information from the ATOM line of the PDB.

class aminoAcid{
public:
  float x[400];
  float y[400];
  float z[400];
  char atomName[1600]; /// 4 times the number of atoms
  char altLoc[400];
  int numAltLoc;
  char possibleLocCode[100];
  char name[4];    // 3 charecter representation of amino acid.
  char charge[3];  // N neutral + positive - negative. -----This field is not used for analysis as it is unreliable.
  int numberOfAtoms; // If number of atoms=0 then residue coordinates not present in PDB file
  int loc;

  aminoAcid(){
    numberOfAtoms=0;
    numAltLoc=0;
    charge[0]='N'; charge[1]='A'; charge[2]='\0'; 
    name[0]='-'; name[1]='-'; name[2]='-'; name[3]='\0';
  }

  ~aminoAcid(){
  }

};

// proteinSequence contains the amino acid sequence of protein. It contains a field resInfo, that 
// contains an array of aminoAcids that holds atom coordinate information.  
// 

class proteinSequence{
public:
  char *sequence;
  int length;
  char seqID;
  aminoAcid *resInfo;
  int currentLoc;

  proteinSequence(){
    sequence=NULL;
    resInfo=NULL;
    currentLoc=0;
  }

  void printSequence(){
    int iLoop;
    cout <<"Length="<<this->length<<" chainID="<<this->seqID<<"\n";
    for(iLoop=0;iLoop<this->length;iLoop++){
      cout <<this->sequence[iLoop*3]<<this->sequence[iLoop*3+1]<<this->sequence[iLoop*3+2]<<" ";
    }
    cout <<"\n";
  }

  void printAtomInfo(){
    int iLoop,jLoop;
    aminoAcid *thsAAList;

    thsAAList=this->resInfo;
    cout <<"Length="<<this->length<<" chainID="<<this->seqID<<"\n";
    for(iLoop=0;iLoop<this->length+10;iLoop++){
      if(thsAAList[iLoop].numberOfAtoms!=0){
	for(jLoop=0;jLoop<thsAAList[iLoop].numberOfAtoms;jLoop++){
	  cout <<thsAAList[iLoop].name<<" "<<thsAAList[iLoop].loc<<" "
	       <<thsAAList[iLoop].atomName[jLoop*4]<<thsAAList[iLoop].atomName[jLoop*4+1]<<thsAAList[iLoop].atomName[jLoop*4+2]
	       <<" "<<thsAAList[iLoop].x[jLoop]<<" "<<thsAAList[iLoop].y[jLoop]<<" "<<thsAAList[iLoop].z[jLoop]<<"\n";
	}
      }
    }
    cout <<"\n";
  }

  ~proteinSequence(){
    if(this->sequence!=NULL)
      delete [] this->sequence;
    if(this->resInfo!=NULL)
      delete [] this->resInfo;
  }
};

// The following class stores HETATM information

class heteroAtom{
public:
  float x;
  float y;
  float z;
  char  atomName[5];
  char name[4];    // 3 charecter representation of amino acid.
  char charge[3];  // N neutral + positive - negative. -----This field is not used for analysis as it is unreliable.
  int index;  // stores the atom index number for this hetatm

  heteroAtom(){
    atomName[0]='-'; atomName[1]='-'; atomName[2]='-'; atomName[3]='-'; atomName[4]='\0';
    charge[0]='N'; charge[1]='A'; charge[2]='\0'; 
    name[0]='-'; name[1]='-'; name[2]='-'; name[3]='\0';
    index=-10;
  }

  ~heteroAtom(){
  }

};


class resultRecord{
public:
  char res1[4];
  char res2[4];
  float distance;
  float angle;
  float angle1;
  float angleP;
  char code[3];
  int loc1;
  int loc2;  
  char chain1;
  char chain2;
  float g_x1;
  float g_y1;
  float g_z1;
  float g_x2;
  float g_y2;
  float g_z2;

  resultRecord(){
    res1[0]='-'; res1[1]='-'; res1[2]='-'; res1[3]='\0';
    res2[0]='-'; res2[1]='-'; res2[2]='-'; res2[3]='\0';
    distance=-100.0;
    angle=1000.0;
    angle1=1000.0;
    angleP=1000.0;
    code[1]='\0'; code[2]='\0'; 
    chain1 = 'z';
    chain2 = 'z';
    g_x1 = -1.000;
    g_y1 = -1.000;
    g_z1 = -1.000;
    g_x2 = -1.000;
    g_y2 = -1.000;
    g_z2 = -1.000;
  }
  
  void setRecord(char *r1,char *r2,float dist,float ang,float ang1,float angP,int l1,int l2,char *cd, char ch1, char ch2, float gx1, float gy1, float gz1, float gx2, float gy2, float gz2){
    memcpy(this->res1,r1,3);
    memcpy(this->res2,r2,3);
    this->distance=dist;
    this->angle=ang;
    this->angle1=ang1;
    this->angleP=angP;
    this->loc1=l1;
    this->loc2=l2;
    this->code[0]=*(cd); this->code[1]=*(cd+1);
    this->chain1 = ch1;
    this->chain2 = ch2;
    this->g_x1 = gx1;
    this->g_y1 = gy1;
    this->g_z1 = gz1;
    this->g_x2 = gx2;
    this->g_y2 = gy2;
    this->g_z2 = gz2;
  }

};


//------------------------------------------------------------------------------------------------------

// the following globals are used to store charge centers that are used for measuring distances
float g_x1[MAX_ALT_LOC],g_y1[MAX_ALT_LOC],g_z1[MAX_ALT_LOC],g_x2[MAX_ALT_LOC],g_y2[MAX_ALT_LOC],g_z2[MAX_ALT_LOC];
int gLocCt;
// the following line has the coordinate of atom in the plane of aromatic ring
float e_x1[MAX_ALT_LOC],e_y1[MAX_ALT_LOC],e_z1[MAX_ALT_LOC];  // this value is set by setPlaneAtomCoord 
int eLocCt;
float h_x1,h_y1,h_z1;  // for hetero atom coordinates
coord P1[MAX_ALT_LOC],P2[MAX_ALT_LOC],P3[MAX_ALT_LOC]; // this value is set by set3PointPlaneCoord. It is used to calculate equation of plane.
coord Q1[MAX_ALT_LOC],Q2[MAX_ALT_LOC],Q3[MAX_ALT_LOC]; // this value is set by set3PointCarboxylPlaneCoord. It is used to calculate equation of plane.
coord CG[MAX_ALT_LOC]; // this value is set by set3PointPlaneCoord. It is used to calculate equation of plane.
int errorCode1[MAX_ALT_LOC],errorCode2[MAX_ALT_LOC],errorCode3[MAX_ALT_LOC],errorCode4[MAX_ALT_LOC]; // this value is set by set3PointPlaneCoord. It is used to calculate equation of plane.


int calculateCenter1(aminoAcid *thsAcid,char *aa,char locCode,int locNumber){
  int iLoop,atomCt=0,eCode=0;
  char *atomType;
  float x[6],y[6],z[6];
  
  g_x1[locNumber]=0.0; 
  g_y1[locNumber]=0.0;
  g_z1[locNumber]=0.0;

  atomType=thsAcid->atomName;
  if((strncmp(aa,"PHE",3)==0)||(strncmp(aa,"TYR",3)==0)){
    for(iLoop=0;iLoop<thsAcid->numberOfAtoms;iLoop++){
      if(((strncmp((atomType+4*iLoop),"CG",2)==0)||(strncmp((atomType+4*iLoop),"CZ",2)==0)||
	  (strncmp((atomType+4*iLoop),"CD1",3)==0)||(strncmp((atomType+4*iLoop),"CD2",3)==0)||
	  (strncmp((atomType+4*iLoop),"CE1",3)==0)||(strncmp((atomType+4*iLoop),"CE2",3)==0))&&
	 (thsAcid->altLoc[iLoop]==locCode)){
	
	x[atomCt]=thsAcid->x[iLoop];
	y[atomCt]=thsAcid->y[iLoop];
	z[atomCt]=thsAcid->z[iLoop];
	atomCt++;
      }
    }
    
    if(atomCt==6){
      for(iLoop=0;iLoop<atomCt;iLoop++){
	g_x1[locNumber]+=x[iLoop];
	g_y1[locNumber]+=y[iLoop];
	g_z1[locNumber]+=z[iLoop];
      }
      g_x1[locNumber]=g_x1[locNumber]/atomCt;
      g_y1[locNumber]=g_y1[locNumber]/atomCt;
      g_z1[locNumber]=g_z1[locNumber]/atomCt;
      eCode=1;
    }else{
      cout <<"Warning: could not locate all the 6 C atoms in benzene ring at residue location "<<thsAcid->loc<<"\n";
      eCode=0;
    }
  }
  return eCode;
}

int calculateTrpCenter1(aminoAcid *thsAcid,char *aa,char locCode,int locNumber){
  int iLoop,atomCt=0,eCode=0;
  char *atomType;
  float x[16],y[16],z[16],wt,wtSum;
  char atomT[16];
  
  g_x1[locNumber]=0.0;
  g_y1[locNumber]=0.0;
  g_z1[locNumber]=0.0;

  atomType=thsAcid->atomName;
  if(strncmp(aa,"TRP",3)==0){
    for(iLoop=0;iLoop<thsAcid->numberOfAtoms;iLoop++){
      if(((strncmp((atomType+4*iLoop),"CG",2)==0)||(strncmp((atomType+4*iLoop),"CH2",3)==0)||
	  (strncmp((atomType+4*iLoop),"CD1",3)==0)||(strncmp((atomType+4*iLoop),"CD2",3)==0)||
	  (strncmp((atomType+4*iLoop),"NE1",3)==0)||(strncmp((atomType+4*iLoop),"CE2",3)==0)||
	  (strncmp((atomType+4*iLoop),"CE3",3)==0)||
	  (strncmp((atomType+4*iLoop),"CZ2",3)==0)||(strncmp((atomType+4*iLoop),"CZ3",3)==0))&&
	 (thsAcid->altLoc[iLoop]==locCode)){	
	x[atomCt]=thsAcid->x[iLoop];
	y[atomCt]=thsAcid->y[iLoop];
	z[atomCt]=thsAcid->z[iLoop];
	atomT[atomCt]=*(atomType+4*iLoop);
	atomCt++;
      }
    }
    wtSum=0.0;
    if(atomCt==9){
      for(iLoop=0;iLoop<atomCt;iLoop++){
	if(atomT[iLoop]=='C'){
	  wt=massC;
	}else if(atomT[iLoop]=='N'){
	  wt=massN;
	}else{
	  eCode=10;
	}
	g_x1[locNumber]+=(x[iLoop])*wt;
	g_y1[locNumber]+=(y[iLoop])*wt;
	g_z1[locNumber]+=(z[iLoop])*wt;
	wtSum+=wt;
      }
      if(wtSum==0.0){
	eCode=10;
      }else{
	g_x1[locNumber]=g_x1[locNumber]/wtSum;
	g_y1[locNumber]=g_y1[locNumber]/wtSum;
	g_z1[locNumber]=g_z1[locNumber]/wtSum;
      }
      if(eCode!=10)
	eCode=1;
      else
	eCode=0;
    }else{
      cout <<"Warning: could not locate all the 9 atoms in TRP ring at residue location "<<thsAcid->loc<<"\n"<<flush;
      eCode=0;
    }
  }
  return eCode;
}

/// The following function calculates the geometric center for carboxylate groups ASP & GLU

int calculateCenter2(aminoAcid *thsAcid,char *aa,char locCode,int locNumber){
  int iLoop,atomCt=0,eCode=0;
  char *atomType;
  float x[6],y[6],z[6];
   
  g_x2[locNumber]=0.0;
  g_y2[locNumber]=0.0;
  g_z2[locNumber]=0.0;
 

  atomType=thsAcid->atomName;
  if((strncmp(aa,"ASP",3)==0)||(strncmp(aa,"GLU",3)==0)){
    for(iLoop=0;iLoop<thsAcid->numberOfAtoms;iLoop++){
      if(thsAcid->altLoc[iLoop]==locCode){	
	if((((strncmp((atomType+4*iLoop),"CG",2)==0)||
	     (strncmp((atomType+4*iLoop),"OD1",3)==0)||
	     (strncmp((atomType+4*iLoop),"OD2",3)==0))&&(strncmp(aa,"ASP",3)==0))||
	   (((strncmp((atomType+4*iLoop),"CD",2)==0)||
	     (strncmp((atomType+4*iLoop),"OE1",3)==0)||
	     (strncmp((atomType+4*iLoop),"OE2",3)==0))&&(strncmp(aa,"GLU",3)==0))){
	
	  x[atomCt]=thsAcid->x[iLoop];
	  y[atomCt]=thsAcid->y[iLoop];
	  z[atomCt]=thsAcid->z[iLoop];
	  atomCt++;
	}  // check for atom type 
      }  // check for alt loc 
    }
    if(atomCt==3){
      for(iLoop=0;iLoop<atomCt;iLoop++){
	g_x2[locNumber]+=x[iLoop];
	g_y2[locNumber]+=y[iLoop];
	g_z2[locNumber]+=z[iLoop];
      }
      g_x2[locNumber]=g_x2[locNumber]/atomCt;
      g_y2[locNumber]=g_y2[locNumber]/atomCt;
      g_z2[locNumber]=g_z2[locNumber]/atomCt;
      eCode=1;
    }else{
      cout <<"Warning: could not locate all the 3 atoms in ASP/GLU side chain at residue location "<<thsAcid->loc<<"\n"<<flush;
      eCode=0;
    }
  }
  return eCode;
}

/// to find the shortest distance among the 2 carboxyl oxygens

void calculateCenter3(aminoAcid *thsAcid,char *aa,char locCode,int locNumber){
  int iLoop,atomCt=0,eCode=0;
  char *atomType;
  float x[6],y[6],z[6];
   
  g_x2[locNumber]=0.0;
  g_y2[locNumber]=0.0;
  g_z2[locNumber]=0.0;
 

  atomType=thsAcid->atomName;
  if((strncmp(aa,"ASP",3)==0)||(strncmp(aa,"GLU",3)==0)){
    for(iLoop=0;iLoop<thsAcid->numberOfAtoms;iLoop++){
      if(thsAcid->altLoc[iLoop]==locCode){	
	if(((strncmp((atomType+4*iLoop),"OD1",3)==0)&&(strncmp(aa,"ASP",3)==0))||
	   ((strncmp((atomType+4*iLoop),"OE1",3)==0)&&(strncmp(aa,"GLU",3)==0))){
	  g_x2[locNumber]=thsAcid->x[iLoop];
	  g_y2[locNumber]=thsAcid->y[iLoop];
	  g_z2[locNumber]=thsAcid->z[iLoop];
	  errorCode3[locNumber]=1;
	}
	if(((strncmp((atomType+4*iLoop),"OD2",3)==0)&&(strncmp(aa,"ASP",3)==0))||
	   ((strncmp((atomType+4*iLoop),"OE2",3)==0)&&(strncmp(aa,"GLU",3)==0))){
	  g_x2[locNumber+1]=thsAcid->x[iLoop];
	  g_y2[locNumber+1]=thsAcid->y[iLoop];
	  g_z2[locNumber+1]=thsAcid->z[iLoop];
	  errorCode3[locNumber+1]=1;
	}
      }  // check for alt loc 
    }
  }
  if((errorCode3[locNumber+1]!=1)||(errorCode3[locNumber]!=1)){
    cout <<"Warning: could not locate at least 1 of the 2 O in ASP/GLU side chain at residue location "<<thsAcid->loc<<"\n"<<flush;
  }
}

/// The following function calculates the geometric center for carboxylate groups ASP & GLU

int calculateCenterOfCharge2(aminoAcid *thsAcid,char *aa,char locCode,int locNumber){
  int iLoop,atomCt=0,eCode=0;
  char *atomType;
  float x[6],y[6],z[6],wt;
   
  g_x2[locNumber]=0.0;
  g_y2[locNumber]=0.0;
  g_z2[locNumber]=0.0;
 

  atomType=thsAcid->atomName;
  if((strncmp(aa,"ASP",3)==0)||(strncmp(aa,"GLU",3)==0)){
    for(iLoop=0;iLoop<thsAcid->numberOfAtoms;iLoop++){
      if(thsAcid->altLoc[iLoop]==locCode){
	int found=0;
	if(((strncmp((atomType+4*iLoop),"OD1",3)==0)||(strncmp((atomType+4*iLoop),"OD2",3)==0))&&(strncmp(aa,"ASP",3)==0)){
	  found=1;
	  wt=chargeO;
	}
	if((strncmp((atomType+4*iLoop),"CG",2)==0)&&(strncmp(aa,"ASP",3)==0)){
	  found=1;
	  wt=chargeC;  
	}
	if((strncmp((atomType+4*iLoop),"CB",2)==0)&&(strncmp(aa,"ASP",3)==0)){
	  found=1;
	  wt=chargeH;  /// treating as H in formic acid
	}
	if(((strncmp((atomType+4*iLoop),"OE1",3)==0)||(strncmp((atomType+4*iLoop),"OE2",3)==0))&&(strncmp(aa,"GLU",3)==0)){
	  found=1;
	  wt=chargeO;
	}
	if((strncmp((atomType+4*iLoop),"CD",2)==0)&&(strncmp(aa,"GLU",3)==0)){
	  found=1;
	  wt=chargeC;  
	}
	if((strncmp((atomType+4*iLoop),"CG",2)==0)&&(strncmp(aa,"GLU",3)==0)){
	  found=1;
	  wt=chargeH;  /// treating as H in formic acid
	}
	
	if(found==1){
	  x[atomCt]=(thsAcid->x[iLoop])*wt;
	  y[atomCt]=(thsAcid->y[iLoop])*wt;
	  z[atomCt]=(thsAcid->z[iLoop])*wt;
	  atomCt++;
	}  // check for atom type 
      }  // check for alt loc 
    }
    if(atomCt==4){
      for(iLoop=0;iLoop<atomCt;iLoop++){
	g_x2[locNumber]+=x[iLoop];
	g_y2[locNumber]+=y[iLoop];
	g_z2[locNumber]+=z[iLoop];
      }
      g_x2[locNumber]=g_x2[locNumber]/(-1.0);  //sum of charges =-1.0
      g_y2[locNumber]=g_y2[locNumber]/(-1.0);
      g_z2[locNumber]=g_z2[locNumber]/(-1.0);
      eCode=1;
    }else{
      cout <<"Warning: could not locate all the 4 atoms in ASP/GLU side chain at residue location "<<thsAcid->loc<<"\n"<<flush;
      eCode=0;
    }
  }
  return eCode;
}

int calculateLysCenter2(aminoAcid *thsAcid,char *aa,char locCode,int locNumber){
  int iLoop,atomCt=0,eCode=0;
  char *atomType;
  float x[6],y[6],z[6],wt;
   
  g_x2[locNumber]=0.0;
  g_y2[locNumber]=0.0;
  g_z2[locNumber]=0.0;

  atomType=thsAcid->atomName;
  if(strncmp(aa,"LYS",3)==0){
    for(iLoop=0;iLoop<thsAcid->numberOfAtoms;iLoop++){
      if(thsAcid->altLoc[iLoop]==locCode){
	if((strncmp((atomType+4*iLoop),"NZ",2)==0)&&(strncmp(aa,"LYS",3)==0)){
	  x[atomCt]=(thsAcid->x[iLoop]);
	  y[atomCt]=(thsAcid->y[iLoop]);
	  z[atomCt]=(thsAcid->z[iLoop]);
	  atomCt++;
	}  // check for atom type 
      }  // check for alt loc 
    }
    if(atomCt==1){
      g_x2[locNumber]=x[0];
      g_y2[locNumber]=y[0];
      g_z2[locNumber]=z[0];
      eCode=1;
    }else{
      cout <<"Warning: could not locate all the 1 atoms in LYS side chain at residue location "<<thsAcid->loc<<"\n"<<flush;
      eCode=0;
    }
  }
  return eCode;
}

int calculateHetAtomCenter2(heteroAtom *thsHetAtm,char *res2,char *atom2,int atomL){
  int iLoop,atomCt=0,eCode=0;
  char *atomType;
  float x[6],y[6],z[6];
   
  h_x1=thsHetAtm->x;
  h_y1=thsHetAtm->y;
  h_z1=thsHetAtm->z;
  eCode=1;
  return eCode;
}

/// The following function assigns the "CG" (c-gamma) or the " edge " atom in the rings to 
/// global variables (e_x1,e_y1,e_z1)

int setPlaneAtomCoord(aminoAcid *thsAcid,char *aa){
  int iLoop,eCode=0;
  char *atomType;

  atomType=thsAcid->atomName;
  if((strncmp(aa,"PHE",3)==0)||(strncmp(aa,"TYR",3)==0)||(strncmp(aa,"TRP",3)==0)){
    for(iLoop=0;iLoop<thsAcid->numberOfAtoms;iLoop++){
      if(strncmp((atomType+4*iLoop),"CG",2)==0){
	e_x1[0]=thsAcid->x[iLoop];
	e_y1[0]=thsAcid->y[iLoop];
	e_z1[0]=thsAcid->z[iLoop];
	eCode=1;
      }
    }
  }
  return eCode;
}

// The following function sets the 3 points in the plane for phe, tyr and trp.
// The variables in which the 3 points are set are P1,P2 and P3.

int set3PointPlaneCoord(aminoAcid *thsAcid,char *aa,char locCode,int locNumber){
  int iLoop,eCode=0;
  char *atomType;

  atomType=thsAcid->atomName;
  if((strncmp(aa,"PHE",3)==0)||(strncmp(aa,"TYR",3)==0)||(strncmp(aa,"TRP",3)==0)){
    for(iLoop=0;iLoop<thsAcid->numberOfAtoms;iLoop++){
      if(thsAcid->altLoc[iLoop]==locCode){	
	if(strncmp((atomType+4*iLoop),"CG",2)==0){
	  P1[locNumber].x=thsAcid->x[iLoop];
	  P1[locNumber].y=thsAcid->y[iLoop];
	  P1[locNumber].z=thsAcid->z[iLoop];
	  CG[locNumber].x=thsAcid->x[iLoop];
	  CG[locNumber].y=thsAcid->y[iLoop];
	  CG[locNumber].z=thsAcid->z[iLoop];
	  eCode+=1;
	}
	if(strncmp((atomType+4*iLoop),"CD1",3)==0){
	  P2[locNumber].x=thsAcid->x[iLoop];
	  P2[locNumber].y=thsAcid->y[iLoop];
	  P2[locNumber].z=thsAcid->z[iLoop];
	  eCode+=1;
	}
	if(strncmp((atomType+4*iLoop),"CD2",3)==0){
	  P3[locNumber].x=thsAcid->x[iLoop];
	  P3[locNumber].y=thsAcid->y[iLoop];
	  P3[locNumber].z=thsAcid->z[iLoop];
	  eCode+=1;
	}
      }
    }
  }
  if(eCode!=3) 
    eCode=0;  //All three points have to be set for error to be non zero
  else
    eCode=1;
  return eCode;
}

// The following function sets the 3 points in the plane for asp and glu.
// The variables in which the 3 points are set are Q1,Q2 and Q3.

int set3PointCarboxylPlaneCoord(aminoAcid *thsAcid,char *aa,char locCode,int locNumber){
  int iLoop,eCode=0;
  char *atomType;

  atomType=thsAcid->atomName;
  if((strncmp(aa,"ASP",3)==0)||(strncmp(aa,"GLU",3)==0)){
    for(iLoop=0;iLoop<thsAcid->numberOfAtoms;iLoop++){
      if(thsAcid->altLoc[iLoop]==locCode){	
	if((strncmp((atomType+4*iLoop),"OD1",3)==0)&&(strncmp(aa,"ASP",3)==0)){
	  Q1[locNumber].x=thsAcid->x[iLoop];
	  Q1[locNumber].y=thsAcid->y[iLoop];
	  Q1[locNumber].z=thsAcid->z[iLoop];
	  eCode+=1;
	}
	if((strncmp((atomType+4*iLoop),"OD2",3)==0)&&(strncmp(aa,"ASP",3)==0)){
	  Q2[locNumber].x=thsAcid->x[iLoop];
	  Q2[locNumber].y=thsAcid->y[iLoop];
	  Q2[locNumber].z=thsAcid->z[iLoop];
	  eCode+=1;
	}
	if((strncmp((atomType+4*iLoop),"CG",2)==0)&&(strncmp(aa,"ASP",3)==0)){
	  Q3[locNumber].x=thsAcid->x[iLoop];
	  Q3[locNumber].y=thsAcid->y[iLoop];
	  Q3[locNumber].z=thsAcid->z[iLoop];
	  eCode+=1;
	}
	if((strncmp((atomType+4*iLoop),"OE1",3)==0)&&(strncmp(aa,"GLU",3)==0)){
	  Q1[locNumber].x=thsAcid->x[iLoop];
	  Q1[locNumber].y=thsAcid->y[iLoop];
	  Q1[locNumber].z=thsAcid->z[iLoop];
	  eCode+=1;
	}
	if((strncmp((atomType+4*iLoop),"OE2",3)==0)&&(strncmp(aa,"GLU",3)==0)){
	  Q2[locNumber].x=thsAcid->x[iLoop];
	  Q2[locNumber].y=thsAcid->y[iLoop];
	  Q2[locNumber].z=thsAcid->z[iLoop];
	  eCode+=1;
	}
	if((strncmp((atomType+4*iLoop),"CD",2)==0)&&(strncmp(aa,"GLU",3)==0)){
	  Q3[locNumber].x=thsAcid->x[iLoop];
	  Q3[locNumber].y=thsAcid->y[iLoop];
	  Q3[locNumber].z=thsAcid->z[iLoop];
	  eCode+=1;
	}
      }
    }
  }
  if(eCode!=3){
    cout <<"Unabel to find all 3 carboxyl group atoms for plane calculation for"<<thsAcid->loc<<"\n"<<flush;
    eCode=0;  //All three points have to be set for error to be non zero
  }else
    eCode=1;
  return eCode;
}

/// The following function calculates the angle between plane and vector (theta).

void calculateAngle(int index1,int index2,float *ang,float *ang1){
  coord P,S,T,X;
  float checkD1,checkD2,checkD3;

  float det=getPlaneEquation(&(P1[index1]),&(P2[index1]),&(P3[index1]),&P);
  S.setCoord(g_x1[index1],g_y1[index1],g_z1[index1]);
  T.setCoord(g_x2[index2],g_y2[index2],g_z2[index2]);
  *ang=getAngleBetweenPlaneAndLine(&P,&S,&T);
  getPlaneProjectCoordinate(&P,&T,-1.0*det,&X);
  *ang1=findAngle(CG[index1].x,CG[index1].y,CG[index1].z,g_x1[index1],g_y1[index1],g_z1[index1],X.x,X.y,X.z);
  checkD1=findDistance(g_x1[index1],g_y1[index1],g_z1[index1],g_x2[index2],g_y2[index2],g_z2[index2]);
  checkD2=findDistance(X.x,X.y,X.z,g_x2[index2],g_y2[index2],g_z2[index2]);
  checkD3=checkD1*float(sin(double(*ang*(3.1412/180.0))));
}

float calculateAngleBetweenPlanes(int index1,int index2){
  coord P,Q;
  
  float det=getPlaneEquation(&(P1[index1]),&(P2[index1]),&(P3[index1]),&P);
  det=getPlaneEquation(&(Q1[index2]),&(Q2[index2]),&(Q3[index2]),&Q);
  float ang=getAngleBetweenPlaneAndPlane(&P,&Q);
  return ang;
}

void processAminoAcid1(aminoAcid *thsAcid,char *aa){
  int iLoop;
  for(iLoop=0;iLoop<thsAcid->numAltLoc;iLoop++){
    if(strncmp(aa,"TRP",3)==0)
      errorCode1[iLoop]=calculateTrpCenter1(thsAcid,aa,thsAcid->possibleLocCode[iLoop],iLoop);
    else
      errorCode1[iLoop]=calculateCenter1(thsAcid,aa,thsAcid->possibleLocCode[iLoop],iLoop);
    errorCode2[iLoop]=set3PointPlaneCoord(thsAcid,aa,thsAcid->possibleLocCode[iLoop],iLoop);
  }
}

void processAminoAcid2(aminoAcid *thsAcid,char *aa){
  int iLoop;
  for(iLoop=0;iLoop<thsAcid->numAltLoc;iLoop++){
    if(rCode)
      errorCode3[iLoop]=calculateCenterOfCharge2(thsAcid,aa,thsAcid->possibleLocCode[iLoop],iLoop);
    else
      calculateCenter3(thsAcid,aa,thsAcid->possibleLocCode[iLoop],(iLoop)*2);

    if((strncmp(aa,"ASP",3)==0)||(strncmp(aa,"GLU",3)==0))
      errorCode4[iLoop]=set3PointCarboxylPlaneCoord(thsAcid,aa,thsAcid->possibleLocCode[iLoop],iLoop);
    else
      errorCode4[iLoop]=0;
  }

}

/// The following subroutine reports the distance between the closest located residues and corresponding angle

resultRecord *findBestInteraction(aminoAcid *thsAcid1,aminoAcid *thsAcid2,resultRecord *opData){
  int totalLoc1=0,totalLoc2=0,iLoop,jLoop;
  float closestDist=10000000000.0,cAngle,cAngle1,cAngleP,angle,angle1,angleP,distance;
  float cg_x1,cg_y1,cg_z1,cg_x2,cg_y2,cg_z2;

  totalLoc1=thsAcid1->numAltLoc;
  totalLoc2=thsAcid2->numAltLoc;
  if((totalLoc1*totalLoc2)==1)
    opData->code[0]='S';
  else
    opData->code[0]='M';

  for(iLoop=0;iLoop<totalLoc1;iLoop++){
    for(jLoop=0;jLoop<totalLoc2;jLoop++){
      if((errorCode1[iLoop]*errorCode2[iLoop]*errorCode3[jLoop])==1){
	distance=findDistance(g_x1[iLoop],g_y1[iLoop],g_z1[iLoop],g_x2[jLoop],g_y2[jLoop],g_z2[jLoop]);
	calculateAngle(iLoop,jLoop,&angle,&angle1);
	if(errorCode4[jLoop]==1)
	  angleP=calculateAngleBetweenPlanes(iLoop,jLoop);
	else
	  angleP=1000.0;
	if(distance<closestDist){
	  closestDist=distance;
	  cAngle=angle;
	  cAngle1=angle1;
	  cAngleP=angleP;
	  cg_x1 = g_x1[iLoop];
	  cg_y1 = g_y1[iLoop];
	  cg_z1 = g_z1[iLoop];
	  cg_x2 = g_x2[jLoop];
	  cg_y2 = g_y2[jLoop];
	  cg_z2 = g_z2[jLoop];
	}
      }
    }
  }
  opData->angle=cAngle;
  opData->angle1=cAngle1;
  opData->distance=closestDist;
  opData->angleP=cAngleP;
  opData->g_x1=cg_x1;
  opData->g_y1=cg_y1;
  opData->g_z1=cg_z1;
  opData->g_x2=cg_x2;
  opData->g_y2=cg_y2;
  opData->g_z2=cg_z2;
  return opData;
}

/// The following subroutine reports the distance between the closest located residues and corresponding angle
/// It is designed to report the closest located oxygen. 

resultRecord *findBestInteraction1(aminoAcid *thsAcid1,aminoAcid *thsAcid2,resultRecord *opData){
  int totalLoc1=0,totalLoc2=0,iLoop,jLoop;
  float closestDist=10000000000.0,cAngle,cAngle1,cAngleP,angle,angle1,angleP,distance;
  float cg_x1,cg_y1,cg_z1,cg_x2,cg_y2,cg_z2;

  totalLoc1=thsAcid1->numAltLoc;
  totalLoc2=thsAcid2->numAltLoc;
  if((totalLoc1*totalLoc2)==1)
    opData->code[0]='S';
  else
    opData->code[0]='M';

  for(iLoop=0;iLoop<totalLoc1;iLoop++){
    for(jLoop=0;jLoop<(totalLoc2)*2;jLoop++){
      if((errorCode1[iLoop]*errorCode2[iLoop]*errorCode3[jLoop])==1){
	distance=findDistance(g_x1[iLoop],g_y1[iLoop],g_z1[iLoop],g_x2[jLoop],g_y2[jLoop],g_z2[jLoop]);
	calculateAngle(iLoop,jLoop,&angle,&angle1);
	if(errorCode4[int(jLoop/2)]==1)
	  angleP=calculateAngleBetweenPlanes(iLoop,int(jLoop/2));
	else
	  angleP=1000.0;
	if(distance<closestDist){
	  closestDist=distance;
	  cAngle=angle;
	  cAngle1=angle1;
	  cAngleP=angleP;
	  cg_x1 = g_x1[iLoop];
	  cg_y1 = g_y1[iLoop];
	  cg_z1 = g_z1[iLoop];
	  cg_x2 = g_x2[jLoop];
	  cg_y2 = g_y2[jLoop];
	  cg_z2 = g_z2[jLoop];
	}
      }
    }
  }
  opData->angle=cAngle;
  opData->angle1=cAngle1;
  opData->distance=closestDist;
  opData->angleP=cAngleP;
  opData->g_x1=cg_x1;
  opData->g_y1=cg_y1;
  opData->g_z1=cg_z1;
  opData->g_x2=cg_x2;
  opData->g_y2=cg_y2;
  opData->g_z2=cg_z2;
  return opData;
}

//////////////////------------------------------------------------------------------
//////------------------------PDB Line Parsing--------------------------------------

int process_SEQRES_Line(char *seqList,int location,char *thsLine){
  char lenC[100],seqC[1000],*tmpStg;
  int seqLen,ct;

  seqLen=charToInt((thsLine+13),4);  // modify based on PDB format

  memcpy(&seqC[0],(thsLine+19),51);  // modify based on PDB format
  seqC[51]='\0';
  
  ct=0;
  tmpStg=strtok(seqC," \t"); //
  do{
    memcpy((seqList+(location+ct)*3),tmpStg,3);
    *(seqList+(location+ct+1)*3)='\0';
    ct++;
    tmpStg=strtok(NULL," \t");
  }while(tmpStg!=NULL);
  return ct;
}

void process_ATOM_Line(char *thsLine,proteinSequence *pdbInfo){
  aminoAcid *thsAAList;
  int currentLoc,aryPosition,iLoop;
  char nameStg[5],*tmpPtr;

  thsAAList=pdbInfo->resInfo;
  currentLoc=charToInt((thsLine+22),4); // modify based on PDB format
  if((currentLoc>=0)&&(*(thsLine+26)==' ')){
    if(currentLoc!=pdbInfo->currentLoc){
      memcpy(&(thsAAList[currentLoc].name[0]),(thsLine+17),3); // modify based on PDB format
      pdbInfo->currentLoc=currentLoc;
      thsAAList[currentLoc].loc=currentLoc;
    }
    thsAAList[currentLoc].x[thsAAList[currentLoc].numberOfAtoms]=charToFloat((thsLine+30),8); // modify based on PDB format
    thsAAList[currentLoc].y[thsAAList[currentLoc].numberOfAtoms]=charToFloat((thsLine+38),8); // modify based on PDB format
    thsAAList[currentLoc].z[thsAAList[currentLoc].numberOfAtoms]=charToFloat((thsLine+46),8); // modify based on PDB format
    if(*(thsLine+16)==' ')  // modify based on PDB format
      thsAAList[currentLoc].altLoc[thsAAList[currentLoc].numberOfAtoms]='X'; // modify based on PDB format
    else{
      thsAAList[currentLoc].altLoc[thsAAList[currentLoc].numberOfAtoms]=*(thsLine+16); // modify based on PDB format
    }
    int isThere=0;
    for(iLoop=0;iLoop<thsAAList[currentLoc].numAltLoc;iLoop++){
      if(thsAAList[currentLoc].possibleLocCode[iLoop]==thsAAList[currentLoc].altLoc[thsAAList[currentLoc].numberOfAtoms])
	isThere=1;
    }
    if(isThere==0){
      thsAAList[currentLoc].possibleLocCode[thsAAList[currentLoc].numAltLoc]=
	thsAAList[currentLoc].altLoc[thsAAList[currentLoc].numberOfAtoms];
      thsAAList[currentLoc].numAltLoc++;
    
    }
    aryPosition=thsAAList[currentLoc].numberOfAtoms*4;
    memcpy(&nameStg[0],(thsLine+12),4);  // modify based on PDB format
    nameStg[4]='\0';
    tmpPtr=strtok(nameStg," \t\n");
    if(tmpPtr!=NULL){
      memcpy(&(thsAAList[currentLoc].atomName[aryPosition]),(tmpPtr),strlen(tmpPtr));
    }
    thsAAList[currentLoc].numberOfAtoms++;
  }
}

void process_HETATM_Line(char *thsLine,heteroAtom *heteroPtr,int position){
  int currentLoc,aryPosition;
  char nameStg[5],*tmpPtr;

  
  memcpy(&(heteroPtr[position].name[0]),(thsLine+17),3); // modify based on PDB format

  heteroPtr[position].x=charToFloat((thsLine+30),8); // modify based on PDB format
  heteroPtr[position].y=charToFloat((thsLine+38),8); // modify based on PDB format
  heteroPtr[position].z=charToFloat((thsLine+46),8); // modify based on PDB format
  memcpy(&nameStg[0],(thsLine+12),4);  // modify based on PDB format
  nameStg[4]='\0';
  tmpPtr=strtok(nameStg," \t\n");
  if(tmpPtr!=NULL){
    memcpy(&(heteroPtr[position].atomName[0]),(tmpPtr),strlen(tmpPtr));
  }
}

int tempCount1234=0;
int getSEQRES_info(char *PDB_File,proteinSequence *pdbInfo){
  char thsLine[LINE_LENGTH];
  char *tmpStg,chainID='-';
  int iLoop,currLength=0,totalLength,chainCt=0;
  ifstream thsPDB,thsPDB1;
  thsPDB.open(PDB_File, ios::in);
  int specCount=0;
  if (!thsPDB.is_open()){ 
    cout << "Error opening PDB file" << PDB_File<<flush;
    exit (1); 
  }

  while (! thsPDB.eof()){
    thsPDB.getline(thsLine,LINE_LENGTH);
    if(strncmp(&thsLine[0],"SEQRES",6)==0){
      if(chainID!=thsLine[11]){
	if(currLength!=0){
	  pdbInfo[chainCt].sequence=new char[currLength*3+1];
	  memcpy(pdbInfo[chainCt].sequence,&tmpSequence[0],(currLength*3));
	  pdbInfo[chainCt].sequence[(currLength*3)]='\0';
	  pdbInfo[chainCt].length=currLength;
	  pdbInfo[chainCt].seqID=chainID;
	  /////	  pdbInfo[chainCt].printSequence();
	  chainCt++;
	}
	currLength=0;
	chainID=thsLine[11];
	totalLength=charToInt(&thsLine[13],4);
      }
      int thsCt=process_SEQRES_Line(&tmpSequence[0],currLength,&thsLine[0]);
      if(tempCount1234==0){
	cout << thsCt << endl;
      }
      currLength+=thsCt;
    }
  }
  tempCount1234++;
  thsPDB.close();
  pdbInfo[chainCt].sequence=new char[currLength*3+1];
  memcpy(pdbInfo[chainCt].sequence,&tmpSequence[0],(currLength*3));
  pdbInfo[chainCt].sequence[(currLength*3)]='\0';
  pdbInfo[chainCt].length=currLength;
  pdbInfo[chainCt].seqID=chainID;
  chainCt++;

  return chainCt;
}

// This function reads all the lines that has ATOM info. from the PDB file. 
// It calls process_ATOM_Line function to get information. It then stores the
// object pdbInfo.

void getATOM_info(char *PDB_File,proteinSequence *pdbInfo){
  char thsLine[LINE_LENGTH];
  char *tmpStg,chainID='-',currentID='-',allChainIDs[1000];
  int iLoop,currLength=0,totalLength,chainCt=-1,seqSize[1000];
  ifstream thsPDB,thsPDB1;
  thsPDB.open(PDB_File, ios::in);
  if (!thsPDB.is_open()){ 
    cout << "Error opening PDB file" << PDB_File<<flush;
    exit (1); 
  }
  while (! thsPDB.eof()){
    thsPDB.getline(thsLine,LINE_LENGTH);
    if(strncmp(&thsLine[0],"ATOM",4)==0){
      if(currentID!=thsLine[21])  // modify based on PDB format
	chainCt++;
      seqSize[chainCt]=charToInt(&thsLine[22],4); // modify based on PDB format
      currentID=thsLine[21]; // modify based on PDB format
      allChainIDs[chainCt]=thsLine[21]; // modify based on PDB format
    }
  }
  thsPDB.close();
  cout <<"Total chains="<<chainCt<<"\n"<<flush;
  for(iLoop=0;iLoop<=chainCt;iLoop++){
    pdbInfo[iLoop].length=seqSize[iLoop];
    if(pdbInfo[iLoop].seqID!=allChainIDs[iLoop])
      cout <<"Warning: chainIDs not matching\n"<<flush;
    pdbInfo[iLoop].resInfo=new aminoAcid[seqSize[iLoop]+10];
  }

  thsPDB1.open(PDB_File, ios::in);
  chainCt=-1;
  if (!thsPDB1.is_open()){ 
    cout << "Error opening PDB file" << PDB_File<<flush;
    exit (1); 
  }
  while (! thsPDB1.eof()){
    thsPDB1.getline(thsLine,LINE_LENGTH);
    if(strncmp(&thsLine[0],"ATOM",4)==0){
      if(chainID!=thsLine[21]){  // modify based on PDB format
	chainCt++;
	chainID=thsLine[21];  // modify based on PDB format
      }
      process_ATOM_Line(&thsLine[0],&(pdbInfo[chainCt]));
      if(strcmp(PDB_File,"/data/AQ/Jason/PDB/1AF7.pdb")==0)
	{
	  cerr << pdbInfo[chainCt].currentLoc << endl;
	}
    }
  }
  thsPDB1.close();
  cout <<"Done reading atom coordinates\n"<<flush;
}

// This function reads all the lines that has ATOM info. from the PDB file. 
// It calls process_ATOM_Line function to get information. It then stores the
// object pdbInfo.

int getHETATM_info(char *PDB_File,heteroAtom *hetAtmList){
  char thsLine[LINE_LENGTH];
  char *tmpStg;
  int iLoop,currLength=0,totalHetAtm,hetAtmCt=0;
  ifstream thsPDB,thsPDB1;
  thsPDB.open(PDB_File, ios::in);

  if (!thsPDB.is_open()){ 
    cout << "Error opening PDB file" << PDB_File<<flush;
    exit (1); 
  }
  totalHetAtm=0;
  while (! thsPDB.eof()){
    thsPDB.getline(thsLine,LINE_LENGTH);
    if(strncmp(&thsLine[0],"HETATM",6)==0){
      totalHetAtm++;
    }
  }
  thsPDB.close();
  
  if(totalHetAtm>MAX_HET_ATOM){
    cout <<"Total Number of Hetero atoms exceeds allocation limit \n"<<flush;
    exit(2);
  }
  
  thsPDB1.open(PDB_File, ios::in);

  if (!thsPDB1.is_open()){ 
    cout << "Error opening PDB file" << PDB_File<<flush;
    exit (1); 
  }
  hetAtmCt=0;
  while (! thsPDB1.eof()){
    thsPDB1.getline(thsLine,LINE_LENGTH);
    if(strncmp(&thsLine[0],"HETATM",6)==0){
      process_HETATM_Line(&thsLine[0],hetAtmList,hetAtmCt);
      hetAtmCt++;
    }
  }
  thsPDB1.close();
  return hetAtmCt;
}


void writeFirstLine(char *fileName){
  int iLoop;
  ofstream opFile;
  opFile.open(fileName, ios::out);

  opFile <<"lineID\tres1\tres2\tdist\tangle\tangleP\tangle1\tloc1\tloc2\tcode\tpdbID\tchain1\tchain2\n";

  opFile.close();
}

void writeRecord(resultRecord *resultArray,int count,char *fileName,char *pdbName){
  int iLoop;
  ofstream opFile;
  opFile.open(fileName, ios::app);

  for(iLoop=0;iLoop<count;iLoop++){
    opFile <<"S\t"<<resultArray[iLoop].res1<<"\t"<<resultArray[iLoop].res2<<"\t"<<resultArray[iLoop].distance
	   <<"\t"<<resultArray[iLoop].angle<<"\t"<<resultArray[iLoop].angleP<<"\t"<<resultArray[iLoop].angle1<<"\t"
	   <<resultArray[iLoop].loc1<<"\t"<<resultArray[iLoop].loc2<<"\t"
	   <<resultArray[iLoop].code<<"\t"<<pdbName<<"\t"<<resultArray[iLoop].chain1<<"\t"<<resultArray[iLoop].chain2<<"\t"<<resultArray[iLoop].g_x1<<"\t"<<resultArray[iLoop].g_y1<<"\t"<<resultArray[iLoop].g_z1<<"\t"<<resultArray[iLoop].g_x2<<"\t"<<resultArray[iLoop].g_y2<<"\t"<<resultArray[iLoop].g_z2<<"\n";
  }

  opFile.close();
}

void searchChainInfo(proteinSequence *pdbInfo,int chainNo1,int chainNo2,char *res1,char *res2,
		     char *name,char *opFile){
  int iLoop,jLoop,recordCt;
  int errorCode1,errorCode2,errorCode3;
  aminoAcid *thsAAList1,*thsAAList2;
  float distance,angle;
  char code[3];
  resultRecord currentResult,*resultPtr;
  resultRecord *interactData= new resultRecord[MAX_RESULTS];
  char ch1;
  char ch2;

  recordCt=0;
  code[0]='I'; code[1]=' '; code[2]='\0';
  if(chainNo1!=chainNo2){
    code[0]='X';  code[1]=' ';
  }
  thsAAList1=pdbInfo[chainNo1].resInfo;
  thsAAList2=pdbInfo[chainNo2].resInfo;
  cout <<"Length="<<pdbInfo[chainNo1].length<<" chainID="<<pdbInfo[chainNo1].seqID<<"\n"<<flush;
  cout <<"Length="<<pdbInfo[chainNo2].length<<" chainID="<<pdbInfo[chainNo2].seqID<<"\n"<<flush;
  cout <<"res1="<<res1<<" res2="<<res2<<"\n"<<flush;

  ch1 = pdbInfo[chainNo1].seqID;
  ch2 = pdbInfo[chainNo2].seqID;

  for(iLoop=0;iLoop<(pdbInfo[chainNo1].length+10);iLoop++){
    if((thsAAList1[iLoop].numberOfAtoms!=0)&&(strncmp(&(thsAAList1[iLoop].name[0]),res1,3)==0)){
      processAminoAcid1(&(thsAAList1[iLoop]),res1);
      for(jLoop=0;jLoop<(pdbInfo[chainNo2].length+10);jLoop++){
	if((thsAAList2[jLoop].numberOfAtoms!=0)&&(strncmp(&(thsAAList2[jLoop].name[0]),res2,3)==0)){
	  processAminoAcid2(&(thsAAList2[jLoop]),res2);
	  if(rCode)
	    resultPtr=findBestInteraction(&(thsAAList1[iLoop]),&(thsAAList2[jLoop]),&currentResult);
	  else
	    resultPtr=findBestInteraction1(&(thsAAList1[iLoop]),&(thsAAList2[jLoop]),&currentResult);
	  code[1]=currentResult.code[0];
	  if(currentResult.distance<15.0){
	    ///////  Subroutines come here
	    interactData[recordCt].setRecord(res1,res2,currentResult.distance,currentResult.angle,currentResult.angle1,currentResult.angleP,
					     iLoop,jLoop,&code[0], ch1, ch2,currentResult.g_x1,currentResult.g_y1,currentResult.g_z1,currentResult.g_x2,currentResult.g_y2,currentResult.g_z2);
	    recordCt++;
	    if(recordCt>=10000){
	      cout <<"Warning: number of interaction records exceeds "<<MAX_RESULTS<<"\n";
	      recordCt=MAX_RESULTS-1;
	    }
	  }             // distance check
	}               // if res2
      }                 // jLoop  
    }                   // if res1
  }                     // iLoop

  cout <<"Total Number of Records="<<recordCt<<"\n";
  writeRecord(interactData,recordCt,opFile,name);
  delete [] interactData;
}


void searchAllChainInteraction(proteinSequence *pdbInfo,int noOfChains, char *name,char *opFile){
  int iLoop,jLoop;

  for(iLoop=0;iLoop<noOfChains;iLoop++){
    for(jLoop=0;jLoop<noOfChains;jLoop++){
      searchChainInfo(pdbInfo,iLoop,jLoop,"TYR","GLU",name,opFile);
      searchChainInfo(pdbInfo,iLoop,jLoop,"TYR","ASP",name,opFile);
      searchChainInfo(pdbInfo,iLoop,jLoop,"PHE","GLU",name,opFile);
      searchChainInfo(pdbInfo,iLoop,jLoop,"PHE","ASP",name,opFile);
      searchChainInfo(pdbInfo,iLoop,jLoop,"TRP","GLU",name,opFile);
      searchChainInfo(pdbInfo,iLoop,jLoop,"TRP","ASP",name,opFile);
    }
  }
}

void processPdbFile(char *name,char *opFile){
  int iLoop;
  proteinSequence *pdbInfo;
  heteroAtom *hetAtmList;
  pdbInfo= new  proteinSequence[20];

  int numOfChains=getSEQRES_info(name,pdbInfo);

  getATOM_info(name,pdbInfo);
  hetAtmList=new heteroAtom[MAX_HET_ATOM];
  int totalHetAtm=getHETATM_info(name,hetAtmList);
  searchAllChainInteraction(pdbInfo,numOfChains,name,opFile);
  //for(iLoop=0;iLoop<numOfChains;iLoop++){
    ///    searchHetAtomInfo(pdbInfo,iLoop,hetAtmList,totalHetAtm,"TYR","HOH","O",1,name,opFile);
    ///    searchHetAtomInfo(pdbInfo,iLoop,hetAtmList,totalHetAtm,"PHE","HOH","O",1,name,opFile);
    ///    searchHetAtomInfo(pdbInfo,iLoop,hetAtmList,totalHetAtm,"TRP","HOH","O",1,name,opFile);
  //}
  if (hetAtmList!=NULL)
    delete [] hetAtmList;
  if (pdbInfo!=NULL)
    delete [] pdbInfo;
}

int processPdbList(char *thsDir,char *pdbList,char *opFile){
  char thsLine[LINE_LENGTH],thsNumber[3];
  char name[MAX_FILE_NAME];
  int iLoop,pdbCount=0,thsLength,pathLength;
  ifstream thsPdbList;

  writeFirstLine(opFile);

  thsPdbList.open(pdbList, ios::in);
  if (!thsPdbList.is_open()){ 
    cout << "Error opening PDB List file " << pdbList << endl;
    exit (1); 
  }
  pathLength=strlen(thsDir);
  //memcpy(&name[0],thsDir,pathLength);
  //name[pathLength]='\0';

  while (! thsPdbList.eof() ){
    thsPdbList.getline(thsLine,LINE_LENGTH);
    thsLength=strlen(&thsLine[0]);
    if((thsLength!=0)&&(!isalpha(thsLine[thsLength-1]))) thsLength--; //remove ^M
    // if( strcmp(thsLine,"1AQZ.pdb")==0 ||
    // 	strcmp(thsLine,"1AIL.pdb")==0 ) continue;
    cout << thsLine<<"\t"<<pdbCount<<"\n"<<flush;
    if(thsLength!=0){
      sprintf(name, "%s/%s",thsDir,thsLine);
      // memcpy(&name[pathLength],&thsLine[0],thsLength);
      // name[pathLength+thsLength]='\0';
      processPdbFile(&name[0],opFile);
      pdbCount++;
    }
  }// end of file 

  thsPdbList.close();

  return pdbCount;
}

int main(int argc, char** argv){
  
  int pdbCt;

  // Parse the command line arguments and store them in an Options class
  Options opt( argc, argv );

  // This is a temp thing since rCode is used elsewhere in the other functions.
  rCode = opt.rCode;

  // process the PDB list file
  pdbCt = processPdbList(opt.pdbfolder, opt.pdbListfile, opt.opFile);

  // if(rCode=='C')  
  //   pdbCt=processPdbList("C:/Scripts/PDB/","C:/Scripts/PDBList.txt","C:/Scripts/Jason/completeResultsCntrOChrg_oxyC_Feb_09_2010.txt");
  // else
  //   pdbCt=processPdbList("C:/Scripts/PDB/","C:/Scripts/PDBList.txt","C:/Scripts/Jason/completeResultsFi_oxyN.txt");
  //pdbCt=processPdbList("/home/51n/courses/lab1/pdbFiles/","completeList.txt","completeResultsFi.txt");


  cout <<"Total PDB count = "<<pdbCt<<"\n";

  return 0;
}
