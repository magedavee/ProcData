#include <iostream>
#include <vector>
#include <string>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1.h>
#include <TClonesArray.h>
#include <TApplication.h>
#include <Event.hh>
#include <sstream>
#define NX 1
#define NY 2
#define NCELL NX*NY
#define ZOFF 450
#define LAT 144.198 
using namespace std;
float ave(vector<float>*);
//float var(vector<float>*);
int main(int argc, char* argv[])
{
    cout<<argc<<" "<<argv[1]<<endl;
    
    if(argc>1)
    {
	for(int i=1;i<argc;++i)
	{
	    //string input("/home/neutrino/mage/CRY_SIM/output/root/");
	    string input("");
	    input.append(argv[i]);
	    cout<<input<<endl;
	    TFile  *file= new TFile(input.c_str());
	    //string name("/home/neutrino/mage/CRY_SIM/output/proc/proc_");
	    string name("");
	    name.append(argv[i]);
	    name.append(".pro");
	    cout<<name<<endl;
	    TTree *tree=(TTree*) file->Get("PG4");
	    IoniClusterEvent *ion=new IoniClusterEvent();
	    ParticleEvent* prim=new ParticleEvent();
	    SecondaryParticleEvent *sp= new SecondaryParticleEvent();
	    tree->SetBranchAddress("Prim",&prim);
	    tree->SetBranchAddress("ScIoni",&ion);
	    tree->SetBranchAddress("SecParticle",&sp);
	    cout<<"Processing "<<argv[i]<<endl;
	    int entries=tree->GetEntries();
	    //int entries=100;
	    TH2F pmtPos("pmtPos","Photon Pos hit",281,-140,140,281,-140+ZOFF+(NY-1)*LAT,140+ZOFF-(NY-1)*LAT);
	    vector<float> * xPos=new vector<float>();
	    vector<float> * yPos=new vector<float>();
	    vector<float> * zPos=new vector<float>();
	    vector<float> * xPos0=new vector<float>();
	    vector<float> * yPos0=new vector<float>();
	    vector<float> * zPos0=new vector<float>();
	    vector<float> * xPos1=new vector<float>();
	    vector<float> * yPos1=new vector<float>();
	    vector<float> * zPos1=new vector<float>();
	    vector<float> * xMom=new vector<float>();
	    vector<float> * yMom=new vector<float>();
	    vector<float> * zMom=new vector<float>();
	    vector<float> * time=new vector<float>();
	    vector<float> * timeL=new vector<float>();
	    vector<float> * timeR=new vector<float>();
	    vector<float> * energy=new vector<float>();
	    vector<int> * PID=new vector<int>();
	    vector<int> * VOL=new vector<int>();
	    vector<int> * EVT=new vector<int>();
	    vector<int> * right=new vector<int>();
	    vector<int> * left=new vector<int>();
	    vector<TH1F*> * lPulse=new vector<TH1F*>();
	    vector<TH1F*> * rPulse=new vector<TH1F*>();
	    TFile *out=new TFile(name.c_str(),"recreate");
	    //TApplication theApp("App",&argc, argv);
	    for(int i=0;i<entries;++i)
	    {
		float percent=i*10000/entries;
		percent=percent/100.0;
		cout<<percent<<" %\n";
		tree->GetEntry(i);
		float x0,y0,z0,x1,y1,z1;
		
		vector<int> * volList=new vector<int>();
		vector<float> * xList[NCELL];
		vector<float> * yList[NCELL];
		vector<float> * zList[NCELL];
		vector<float> * tLList[NCELL];
		vector<float> * tRList[NCELL];
		for(int j =0;j<NCELL;++j)
		{
		    xList[j]=new vector<float>();
		    yList[j]=new vector<float>();
		    zList[j]=new vector<float>();
		    tLList[j]=new vector<float>();
		    tRList[j]=new vector<float>();
		}
		
		TClonesArray * photon=sp->particles;
		TClonesArray * particle=prim->particles;
		TClonesArray * clusts=ion->clusts;
		
		IoniCluster* v = (IoniCluster*) clusts->At(0);
		float t0=v->t;
		float e=ion->EIoni;	
		int det=sp->nParticlesDet;
		int part=prim->nParticles;
		if(part==1)
		{
		    ParticleVertex* pv=(ParticleVertex*) particle->At(0);
		    int p=pv->PID;
		    float mx=pv->p[0];
		    float my=pv->p[1];
		    float mz=pv->p[2];
		    if(mz<-.12228)
		    {
			string lName("left");
			string rName("right");
			ostringstream oss;
			oss << i;
			lName+=oss.str();
			rName+=oss.str();
			SecondaryParticleVertex* temp=(SecondaryParticleVertex*) photon->At(0);
			float max=temp->t-t0;
			float min=temp->t-t0;
			for(int j=0;j<det;++j)
			{
			    SecondaryParticleVertex* pmt=(SecondaryParticleVertex*) photon->At(j);
			    float t=pmt->t-t0;
			    if(min>t)
			    {
				min=t;
			    }
			    if(max<t)
			    {
				max=t;
			    }
			    float y=pmt->x[1];
			    float x=pmt->x[0];
			    float z=pmt->x[2]+(NY-1)*LAT+ZOFF;
			    int cell=z/LAT;
			    pmtPos.Fill(x,z);
			    if(y<0)
			    {
				    tLList[cell]->push_back(t);
			    }
			    else
			    {
				    tRList[cell]->push_back(t);
			    }
			}
			float tLAve[NCELL];
			float tRAve[NCELL];
			for(int j=0;j<NCELL;++j)
			{
			    tLAve[j]=ave(tLList[j]);
			    tRAve[j]=ave(tRList[j]);
			}
			int photoL[NCELL];
			int photoR[NCELL];
			for(int j=0;j<NCELL;++j)
			{
			    photoL[j]=tLList[j]->size();
			    photoR[j]=tRList[j]->size();
			}
			float dT[NCELL];
			for(int j=0;j<NCELL;++j)
			{
			   dT[j] =tLAve[j]-tRAve[j];
			}
			bool test=true;
			for(int j=0;j<NCELL;++j)
			{
			   if ((dT[j]*dT[j])>(100*100))
			   {
			       test=false;
			   }
			}
			if(test)
			{
			    int num=ion->nIoniClusts;
			    for(int j=0;j<num;++j)
			    {
				    IoniCluster* vert = (IoniCluster*) clusts->At(j);
				    float x=vert->x[0];
				    float y=vert->x[2];
				    float z=vert->x[1];
				    if(j==0)
				    {
					x0=x;
					y0=y;
					z0=z;
				    }
				    x1=x;
				    y1=y;
				    z1=z;
				    int pid=vert->PID;
				    int vol=vert->vol;
				    if(pid==p&& vol>=0)
				    {
					if(vol<NCELL)
					{
					    xList[vol]->push_back(x);
					    yList[vol]->push_back(y);
					    zList[vol]->push_back(z);
					    volList->push_back(vol);
					}
					else
					{
					    cout<<vol<<" is greater than NCELL value "<<NCELL<<" and needs to be changed\n";
					}
				    }
			    }
			    for(int j =0;j<NCELL;++j)
			    {
				if(xList[j]->size()!=0)
				{
				    float xAve=ave(xList[j]);
				    float yAve=ave(yList[j]);
				    float zAve=ave(zList[j]);
				    float dt=dT[j];
				    float tl=tLAve[j];
				    float tr=tRAve[j];
				    int Left=photoL[j];
				    int Right=photoR[j];
				    left->push_back(Left);
				    right->push_back(Right);
				    xPos->push_back(xAve);
				    yPos->push_back(yAve);
				    zPos->push_back(zAve);
				    xPos0->push_back(x0);
				    yPos0->push_back(y0);
				    zPos0->push_back(z0);
				    xPos1->push_back(x1);
				    yPos1->push_back(y1);
				    zPos1->push_back(z1);
				    time->push_back(dt);
				    timeL->push_back(tl);
				    timeR->push_back(tr);
				    energy->push_back(e);
				    PID->push_back(p);
				    VOL->push_back(j);
				    EVT->push_back(i);
				    xMom->push_back(mx);
				    yMom->push_back(my);
				    zMom->push_back(mz);
				    //lPulse->push_back(left);
				    //rPulse->push_back(right);
				}
			    }
			} 
		    }
		    else
			cout<<mz<<endl;
		}
	    }
	    //pmtPos.Draw("colz");
	    //char q;
	    //cin>>q;
	    //theApp.Run();
	    cout<<"creating "<<name<<endl;
	    TTree treeOut("photon_Data","process photon data");
	    float x,y,z0,x0,y0,z1,x1,y1,z,t,e,mX,mY,mZ,tl,tr;
	    int p,v,ev,l,r;
	    treeOut.Branch("x",&x,"x/F");
	    treeOut.Branch("y",&y,"y/F");
	    treeOut.Branch("z",&z,"z/F");
	    treeOut.Branch("x0",&x0,"x0/F");
	    treeOut.Branch("y0",&y0,"y0/F");
	    treeOut.Branch("z0",&z0,"z0/F");
	    treeOut.Branch("x1",&x1,"x1/F");
	    treeOut.Branch("y1",&y1,"y1/F");
	    treeOut.Branch("z1",&z1,"z1/F");
	    treeOut.Branch("t",&t,"t/F");
	    treeOut.Branch("tl",&tl,"tl/F");
	    treeOut.Branch("tr",&tr,"tr/F");
	    treeOut.Branch("e",&e,"e/F");
	    treeOut.Branch("momX",&mX,"mX/F");
	    treeOut.Branch("momY",&mY,"mY/F");
	    treeOut.Branch("momZ",&mZ,"mZ/F");
	    treeOut.Branch("pid",&p,"p/I");
	    treeOut.Branch("vol",&v,"v/I");
	    treeOut.Branch("evt",&ev,"ev/I");
	    treeOut.Branch("left",&l,"l/I");
	    treeOut.Branch("right",&r,"r/I");
	    int evt=xPos->size();
	    for(int j=0;j<evt;++j)
	    {
		cout<<"writing "<<j<<" event out of "<<evt<<endl;
		x=xPos->at(j);
		y=yPos->at(j);
		z=zPos->at(j);
		x0=xPos0->at(j);
		y0=yPos0->at(j);
		z0=zPos0->at(j);
		x1=xPos1->at(j);
		y1=yPos1->at(j);
		z1=zPos1->at(j);
		t=time->at(j);
		tl=timeL->at(j);
		tr=timeR->at(j);
		e=energy->at(j);
		p=PID->at(j);
		v=VOL->at(j);
		ev=EVT->at(j);
		mX=xMom->at(j);
		mY=yMom->at(j);
		mZ=zMom->at(j);
		l=left->at(j);
		r=right->at(j);
		treeOut.Fill();
		//lPulse->at(j)->Write();
		//rPulse->at(j)->Write();
	    }
	    cout<<"writing file\n";
	    treeOut.Write();
	    out->Close();

	    
	}
    }
    else
    {
	cout<<"no input file\n";
    }
    return 0;
}

float ave(vector<float>* x)
{
    float sum=0;
    float size=(float)x->size();
    //cout<<"size :"<<size<<endl;
    if(size!=0)
    {
	for(auto it=x->begin();it!=x->end();++it)
	{
	    float num=*it;
	    if(!isnan(num))
		sum+=num;
	    else
		cout<<"not a number, read ignored "<<endl;
	}
	return sum/size;
    }
    else
	return -1;
}

//float var(vector<float>* x)
//{
    //float sum=0;
    //float size=(float)x->size();
    ////cout<<"size :"<<size<<endl;
    //if(size!=0)
    //{
	//for(auto it=x->begin();it!=x->end();++it)
	//{
	    //float num=*it;
	    //if(!isnan(num))
		//sum+=num;
	    //else
		//cout<<"not a number, read ignored "<<endl;
	//}
	//return sum/size;
    //}
    //else
	//return -1;
//}
