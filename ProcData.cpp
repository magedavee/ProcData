#include <iostream>
#include <vector>
#include <string>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1.h>
#include <TClonesArray.h>
#include <Event.hh>
#include <sstream>
#define NX 14
#define NY 10
#define NCELL NX*NY
using namespace std;
float ave(vector<float>*);
//float var(vector<float>*);
int main(int argc, char* argv[])
{
    
    if(argc>1)
    {
	for(int i=1;i<argc;++i)
	{
	    //string input("/home/neutrino/mage/CRY_SIM/output/root/");
	    string input("");
	    input.append(argv[i]);
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
	    vector<float> * xPos=new vector<float>();
	    vector<float> * yPos=new vector<float>();
	    vector<float> * zPos=new vector<float>();
	    vector<float> * time=new vector<float>();
	    vector<float> * energy=new vector<float>();
	    vector<int> * PID=new vector<int>();
	    vector<int> * VOL=new vector<int>();
	    vector<int> * EVT=new vector<int>();
	    vector<TH1F*> * lPulse=new vector<TH1F*>();
	    vector<TH1F*> * rPulse=new vector<TH1F*>();
	    TFile *out=new TFile(name.c_str(),"recreate");
	    for(int i=0;i<entries;++i)
	    {
		float percent=i*10000/entries;
		percent=percent/100.0;
		cout<<percent<<" %\n";
		tree->GetEntry(i);

		
		vector<int> * volList=new vector<int>();
		vector<float> * xList[NCELL];
		vector<float> * yList[NCELL];
		vector<float> * zList[NCELL];
		for(int j =0;j<NCELL;++j)
		{
		    xList[j]=new vector<float>();
		    yList[j]=new vector<float>();
		    zList[j]=new vector<float>();
		}
		//vector<float> * xList=new vector<float>();
		//vector<float> * yList=new vector<float>();
		//vector<float> * zList=new vector<float>();
		vector<float> * tLList=new vector<float>();
		vector<float> * tRList=new vector<float>();
		
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
		    string lName("left");
		    string rName("right");
		    ostringstream oss;
		    oss << i;
		    lName+=oss.str();
		    rName+=oss.str();
		    //cout<<lName<<endl;
		    SecondaryParticleVertex* temp=(SecondaryParticleVertex*) photon->At(0);
		    float max=temp->t-t0;
		    float min=temp->t-t0;
		    //cout<<"max "<<temp->t<<endl;
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
			//cout<<t<<endl;
			float y=pmt->x[1];
			if(y<0)
			{
				tLList->push_back(t);
			}
			else
			{
				tRList->push_back(t);
			}
		    }
		    //int Min=(int)min;
		    //int Max=(int)max;
		    //cout<<"min "<<Min;
		    //TH1F* left= new TH1F(lName.c_str(),lName.c_str(),Max-Min,Min,Max);
		    //TH1F* right= new TH1F(lName.c_str(),lName.c_str(),Max-Min,Min,Max);
		    //for(int j=0;j<det;++j)
		    //{
			//SecondaryParticleVertex* pmt=(SecondaryParticleVertex*) photon->At(j);
			//float t=pmt->t-t0;
			//float y=pmt->x[1];
			//if(y<0)
				//left->Fill(t);
			//else
				//right->Fill(t);
		    //}
		    float tLAve=ave(tLList);
		    float tRAve=ave(tRList);
		    float dT=tLAve-tRAve;
		    if((dT*dT)<(100*100))
		    {
			int num=ion->nIoniClusts;
			for(int j=0;j<num;++j)
			{
				IoniCluster* vert = (IoniCluster*) clusts->At(j);
				float x=vert->x[0];
				float y=vert->x[2];
				float z=vert->x[1];
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
				//cout<<"x pos:"<<xPos<<endl;
				//cout<<"y pos:"<<yPos<<endl;
				//cout<<"z pos:"<<zPos<<endl;
				xPos->push_back(xAve);
				yPos->push_back(yAve);
				zPos->push_back(zAve);
				time->push_back(dT);
				energy->push_back(e);
				PID->push_back(p);
				VOL->push_back(j);
				EVT->push_back(i);
				//lPulse->push_back(left);
				//rPulse->push_back(right);
			    }
			}
		    
		    }
		}
	    }
	    
	    cout<<"creating "<<name<<endl;
	    TTree treeOut("photon_Data","process photon data");
	    float x,y,z,t,e;
	    int p,v,ev;
	    treeOut.Branch("x",&x,"x/F");
	    treeOut.Branch("y",&y,"y/F");
	    treeOut.Branch("z",&z,"z/F");
	    treeOut.Branch("t",&t,"t/F");
	    treeOut.Branch("e",&e,"e/F");
	    treeOut.Branch("pid",&p,"p/I");
	    treeOut.Branch("vol",&v,"v/I");
	    treeOut.Branch("evt",&ev,"ev/I");
	    int evt=xPos->size();
	    for(int j=0;j<evt;++j)
	    {
		cout<<"writing "<<j<<" event out of "<<evt<<endl;
		x=xPos->at(j);
		y=yPos->at(j);
		z=zPos->at(j);
		t=time->at(j);
		e=energy->at(j);
		p=PID->at(j);
		v=VOL->at(j);
		ev=EVT->at(j);
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
