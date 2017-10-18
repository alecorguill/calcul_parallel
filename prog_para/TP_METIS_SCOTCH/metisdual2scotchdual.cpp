#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <list>

using namespace std;


// Un élément est une liste d'entiers (les voisins)
typedef list<int> Element;

// Un maillage dual contient une liste d'éléments
class MaillageDual
{
	public:

	MaillageDual()
	{
	}
	~MaillageDual()
	{
	}
	
	list<Element> _elements;
	int _n2;
	
	// Remplir le maillage dual depuis un fichier
	int lire( string fileName )
	{
		ifstream fichier(fileName.data(), ios::in);
	
		if( !fichier.good() )
		{
			cerr << "Ne peut pas ouvrir " << fileName << " pour lecture !";
			return -1;
		}
	
		cout << "Lecture du fichier " << fileName << "... ";
		
		string ligne;

		int nelements;
		
		getline(fichier,ligne);
		if( sscanf( ligne.data(), "%d %d", &nelements, &this->_n2 ) != 2 )
		{
			cerr << "Ne trouve pas le nombre d'elements !";
			return -2;
		}
		cout << nelements << " elements... ";

		// Traite chaque élément
		string delim = " ";
		for( int i=0 ; i<nelements ; i++ )
		{
			getline(fichier,ligne);
			
			Element e;
			
			int start = 0;
			int end = ligne.find(delim);
			
			while( start < ligne.length() )
			{
				int n; // numéro lu
				int res; // nb int lus
				
				res = sscanf( ligne.substr(start, end-start ).data(), "%d", &n );
				
				if( res == 1 )
				{
					e.push_back( n );
				}
							
				if( end != -1 )
				{
					start = end + delim.length();
					end = ligne.find(delim, start);
				}
				else
				{
					start = ligne.length()+1;
				}
				
			}
			this->_elements.push_back(e);
			
		}
		
		cout << "lecture terminee !\n";
		
		fichier.close();
		
		return 0;
	}
	
	// Remplir le maillage dual au format Scotch dans un fichier
	int ecrireScotch( string fileName )
	{
	
		ofstream fichier(fileName.data(), ios::out);
	
		if( !fichier.good() )
		{
			cerr << "Ne peut pas ouvrir " << fileName << " pour écriture !";
			return -1;
		}
	
		cout << "Ecriture du fichier " << fileName << "... ";
		
		string ligne;

		fichier << "0\n";
		fichier << this->_elements.size() << " " << this->_n2*2 << "\n";
		fichier << "0 000\n";
		
		list<Element>::iterator it;
		for( it=this->_elements.begin() ; it!=this->_elements.end() ; it++ )
		{
			fichier << it->size();
			
			Element::iterator ite;
			for( ite=it->begin(); ite!=it->end() ; ite++ )
			{
				fichier << " " << (*ite)-1;
			}
			fichier << "\n";
		}

		cout << "ecriture terminee !\n";
		
		fichier.close();
		
		return 0;
	}
	
	
};

int main( int argc, char **argv )
{
	
	if( argc != 3 )
	{
		std::cerr << "usage: " << argv[0] << " <metis dual mesh file> <scotch dual mesh file>";
		return -1;
	}
	
	string inputFileName = argv[1];
	string outputFileName = argv[2];
	
	MaillageDual maillage;
	
	maillage.lire( inputFileName );
	maillage.ecrireScotch( outputFileName );

		return 0;
}

