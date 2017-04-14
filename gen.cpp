#include <fstream>

int main()
{
  std::ofstream out;
  out.open("test4");

  for(int i=1; i<=50; ++i)
  {
    for(int j=1; j<=9; ++j)
    {
      out << (rand()%150) << ",";
    }
    out<<(rand()%150)<<"\n";
  }
  
  out.close();
}
