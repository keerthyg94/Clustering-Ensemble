import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;


public class ensemble {
	
	public static void main(String args[]) throws Exception
	{  
		String filename=null;
		String op_dc=null,op_hc=null,op_pc=null;
		
		Scanner scan= new Scanner(System.in);
		
		int mGenes=0;
		System.out.println(" Clustering Ensemble...");
		System.out.println(" Enter the dataset..");
		BufferedReader reader=new BufferedReader(new InputStreamReader(System.in));
        filename = reader.readLine();
		
		BufferedReader readfile = new BufferedReader(new FileReader(filename));
		while(readfile.readLine()!=null)	
		{
			mGenes++;
		}
		BufferedReader readfile1 = new BufferedReader(new FileReader(filename));
		String read;
		read=readfile1.readLine();
		String tokens[]=read.split("\\t");
		int columns =tokens.length-2;
		//System.out.println("Total no of transactions ="+mGenes);
		int groundtruth[] = new int[mGenes];
		int value_dc[] = new int[mGenes];
		int value_pc[] = new int[mGenes];
		int value_hc[] = new int[mGenes];
		double data_matrix[][] = new double[mGenes][columns];
		BufferedReader readfile2 = new BufferedReader(new FileReader(filename));
		String read1;
		int row1=0;
		int column1=0;
		while ((read1 = readfile2.readLine()) != null) 
		{
			//read1=readfile2.readLine();
			String tokens1[]=read1.split("\\t");
			column1=0;
			for(int i=0;i<=1;i++)
			{
				String geneid=tokens1[0];
				String truth=tokens1[1];
				int gid=Integer.parseInt(geneid);
				int gtruth=Integer.parseInt(truth);
				groundtruth[row1]= gtruth;
			}
			for(int i=2;i<tokens1.length;i++)
			{
				double value=Double.parseDouble(tokens1[i]);
				//System.out.println(value);
				data_matrix[row1][column1]=value;
				column1++;
				
			}
			row1++;
		}
		
		
		double distance_matrix[][] = new double[mGenes][mGenes];
		//double new_distance_matrix[][]= new double[mGenes][mGenes];
		
		System.out.println();
		for(int i=0;i<row1;i++)
		{
			for(int j=0;j<row1;j++)
			{
				if(i==j)
				{
					distance_matrix[i][j]=0.0;
					//new_distance_matrix[i][j]=0.0;
				}
				else
				{	double sum=0;
				    double temp=0;
					for(int k=0;k<column1;k++)
					{  
						temp=data_matrix[i][k]-data_matrix[j][k];
						double sq_temp=temp*temp;
						sum=sum+sq_temp;
					}
					double root_temp=Math.sqrt(sum);
					distance_matrix[i][j]=root_temp;
					//new_distance_matrix[i][j] = root_temp;
					
				}
			}
			//System.out.println();
		}//end of for
		
		
		//System.out.println(" Enter the o/p file of density based clustering");
		if(filename.equals("cho.txt"))
		{
		op_dc="cho_dc.txt";
		BufferedReader readfile_dc = new BufferedReader(new FileReader(op_dc));
		String read_dc;
		int row_dc=0;
		while ((read_dc = readfile_dc.readLine()) != null) 
		{
			//read1=readfile2.readLine();
			String tokens1[]=read_dc.split("\\t");
			String truth=tokens1[1];
			int gtruth=Integer.parseInt(truth);
			value_dc[row_dc]= gtruth;
			row_dc++;
		}
		
		//System.out.println(" Enter the o/p file of hierarchical clustering");
		op_hc="cho_hc.txt";
		BufferedReader readfile_hc = new BufferedReader(new FileReader(op_hc));
		String read_hc;
		int row_hc=0;
		while ((read_hc = readfile_hc.readLine()) != null) 
		{
			//read1=readfile2.readLine();
			String tokens1[]=read_hc.split("\\t");
			String truth=tokens1[1];
			int gtruth=Integer.parseInt(truth);
			value_hc[row_hc]= gtruth;
			row_hc++;
		}
		
		//System.out.println(" Enter the o/p file of partitional clustering");
		op_pc="cho_pc.txt";
		BufferedReader readfile_pc = new BufferedReader(new FileReader(op_pc));
		String read_pc;
		int row_pc=0;
		while ((read_pc = readfile_pc.readLine()) != null) 
		{
			//read1=readfile2.readLine();
			String tokens1[]=read_pc.split("\\t");
			String truth=tokens1[1];
			int gtruth=Integer.parseInt(truth);
			value_pc[row_pc]= gtruth;
			row_pc++;
		}
	}
		else if(filename.equals("iyer.txt"))
		{
			op_dc="iyer_dc.txt";
			BufferedReader readfile_dc = new BufferedReader(new FileReader(op_dc));
			String read_dc;
			int row_dc=0;
			while ((read_dc = readfile_dc.readLine()) != null) 
			{
				//read1=readfile2.readLine();
				String tokens1[]=read_dc.split("\\t");
				String truth=tokens1[1];
				int gtruth=Integer.parseInt(truth);
				value_dc[row_dc]= gtruth;
				row_dc++;
			}
			
			//System.out.println(" Enter the o/p file of hierarchical clustering");
			op_hc="iyer_hc.txt";
			BufferedReader readfile_hc = new BufferedReader(new FileReader(op_hc));
			String read_hc;
			int row_hc=0;
			while ((read_hc = readfile_hc.readLine()) != null) 
			{
				//read1=readfile2.readLine();
				String tokens1[]=read_hc.split("\\t");
				String truth=tokens1[1];
				int gtruth=Integer.parseInt(truth);
				value_hc[row_hc]= gtruth;
				row_hc++;
			}
			
			//System.out.println(" Enter the o/p file of partitional clustering");
			op_pc="iyer_pc.txt";
			BufferedReader readfile_pc = new BufferedReader(new FileReader(op_pc));
			String read_pc;
			int row_pc=0;
			while ((read_pc = readfile_pc.readLine()) != null) 
			{
				//read1=readfile2.readLine();
				String tokens1[]=read_pc.split("\\t");
				String truth=tokens1[1];
				int gtruth=Integer.parseInt(truth);
				value_pc[row_pc]= gtruth;
				row_pc++;
			}
		}
		else
		{
			System.out.println(" Incorrect file entered..");
			System.exit(0);
		}
		int groundvalue = num_elements(groundtruth);
		System.out.println(" Ground value = "+groundvalue);
		
		relable_hierrarchichal(value_hc, mGenes, groundvalue);
		/*System.out.println(" The output of the file is as follows");
		for(int i=0;i<mGenes;i++)
		{
			System.out.println(i + "   "+ value_hc[i]);
		}*/
		
		
		
		
		int return_hc[]= new int[groundvalue];
		return_hc=valueArray(value_hc, mGenes, groundvalue);
		for(int i=0;i<return_hc.length;i++)
		{
			System.out.println(i + "   "+return_hc[i]);
		}
		
		int return_pc[] = new int[groundvalue];
		return_pc=valueArray(value_pc, mGenes, groundvalue);
		for(int i=0;i<return_pc.length;i++)
		{
			System.out.println(i + "   "+return_pc[i]);
		}
		int radius=groundvalue;
		int totalvalue= Factorial(groundvalue);
		int numTotal=totalvalue;
		int b[]=new int[groundvalue];
		int indices[];
		for(int i=0;i<b.length;i++)
		{
			b[i]=i+1;
		}
		int cluster_ensemble[] = new int[mGenes];
		int foundtruth[] = new int[mGenes];
		for(int i=0;i<mGenes;i++)
		{
			if(value_pc[i]==value_hc[i])
			{
				cluster_ensemble[i] = value_hc[i];
			}
			else if(value_hc[i]==value_dc[i])
			{
				cluster_ensemble[i]= value_hc[i];
			}
			else if(value_pc[i]==value_dc[i])
			{
				cluster_ensemble[i] = value_pc[i];
			}
			else
			{
				cluster_ensemble[i] = value_pc[i];
			}
		}
		
		
		/*while(totalvalue>0)
		{   
			indices = DifferentCombinations(totalvalue,radius,groundvalue,numTotal,b);
			for(int l=0;l<indices.length;l++)
			{
				System.out.print(indices[l] + "  ");
			}
			System.out.println();
			
			b=indices;
			totalvalue--;
		}*/
		
		System.out.println(" The clustering ensembling result is as follows..");
		for(int i=0;i<mGenes;i++)
		{  foundtruth[i] =cluster_ensemble[i];
			System.out.println("Gene Id="+(i +1)+ " Cluster number= " + cluster_ensemble[i]);
		}
		
		
		//Incident_matrix calculation
				int incident_matrixP[][] = new int[mGenes][mGenes];
				int incident_matrixC[][] = new int[mGenes][mGenes];
				int ss=0,dd=0,sd=0,ds=0;
				for(int i=0;i<mGenes;i++)
				{
					for(int j=0; j<mGenes; j++)
					{
						if(groundtruth[i]==groundtruth[j])
						{
							incident_matrixP[i][j] =1;
						}
						
						else
						{
							incident_matrixP[i][j]=0;
						}
						
						if(foundtruth[i]==foundtruth[j])
						{
							incident_matrixC[i][j]=1;
						}
						else
						{
							incident_matrixC[i][j]=0;
						}
						
					}
				}//end of for
				
			for(int i=0;i<mGenes;i++)	
			{
				for(int j=0;j<mGenes;j++)
				{
					if(incident_matrixP[i][j]==incident_matrixC[i][j])
					{
						if(incident_matrixP[i][j]==1)
						{
							ss++;
						}
						else
						{
							dd++;
						}
					}
					else
					{
						if(incident_matrixP[i][j]==1)
						{
							ds++;
						}
						else
						{
							sd++;
						}
					}
				}
			}
				
			System.out.println("SS = "+ss);
			System.out.println("DD = "+dd);
			System.out.println("SD = "+sd);
			System.out.println("DS = "+ds);
			
			double rand_index;
			double jaccard_coef;
			double ss1=(double)ss;
			double dd1 = (double)dd;
			double sd1= (double)sd;
			double ds1 = (double)ds;
			
			rand_index =(ss1+dd1)/(ss1+dd1+sd1+ds1);
			jaccard_coef=(ss1)/(ss1+sd1+ds1);
			
			System.out.println(" Rand index = "+rand_index);
			System.out.println(" Jacard Coefficient = "+jaccard_coef);
			
			//internal_index calculation
			double cincident_matrix[][] = new double[mGenes][mGenes];
			for(int i=0;i<mGenes;i++)
			{
				for(int j=0; j<mGenes;j++)
				{
					cincident_matrix[i][j] =(double)incident_matrixC[i][j];
				}
			}
			
			double dmean = mean(distance_matrix,mGenes);
			double cmean = mean(cincident_matrix,mGenes);
			double nvar=numerical_variance(distance_matrix, cincident_matrix, mGenes);
			//System.out.println(" Nvar = "+nvar);
			double dvar=Math.sqrt(variance(distance_matrix,mGenes));
			//System.out.println(" Dvar = "+dvar);
			double cvar = Math.sqrt(variance(cincident_matrix, mGenes));
			//System.out.println(" Cvar = "+cvar);
			double correlation = Math.abs(nvar/(dvar*cvar));
			System.out.println(" Correlation of incident matrix and distance matrix = "+correlation);
		
	}
	
	public static void relable_hierrarchichal(int arr[], int num,int val)
	{  
		Set<Integer> newset = new HashSet<Integer>();
		for(int i=1;i<=val;i++)
		{
			newset.add(i);
		}
		int clusternum=2;
		for(int i=0;i<num;i++)
		{
			if(!newset.contains(arr[i]))
			{  //System.out.println("arr[i]="+arr[i]);
				for(int j=i+1;j<num;j++)
				{
					if(arr[j]==arr[i])
					{  
						arr[j]=clusternum;
						//System.out.println("arr[j]= "+arr[j]);
					}
				}
				arr[i]=clusternum;
				clusternum++;
				System.out.println(" cluster number = "+clusternum);
			}
		}
			
	}
	
	public static int vote(int arr1[],int arr2[], int arr3[],int num)
	{
		int count =0;
		for(int i=0;i<num;i++)
		{
			
			if((arr1[i]==arr2[i])&&(arr1[i]==arr3[i]))
			{
				count++;
			}
		}
		return count;
	}
	
	
	public static int[] valueArray(int arr[],int num, int val)
	{
		Set<Integer> newset = new HashSet<Integer>();
		Set<Integer> newset1 = new HashSet<Integer>();
		newset1.add(-1);
		int count=0;
		int returnarr[] = new int[val];
		for(int i=1;i<=val;i++)
		{
			newset.add(i);
		}
		for(int i=0;i<num;i++)
		{
			if(!newset1.contains(arr[i]))
			{
				returnarr[count]=i;
				count++;
				newset1.add(arr[i]);
			}
		}
		
		return returnarr;
	}
	
	
	//function for different elements in groundtruth
	public static int num_elements(int[] arr)
	{
		Set<Integer> newset = new HashSet<Integer>();
		for(int element: arr)
		{
			newset.add(element);
		}
		if(newset.contains(-1))
		{
			return newset.size()-1;
		}
		else
		{
			return newset.size();
		}
		
	}
	
	public static int Factorial(int n)
	{
		int  fact=1;
		for (int i = n; i > 1; i--) 
		{
			fact = fact*i;
		}
		return fact;
	}
	
	public static int[] DifferentCombinations(int total,int radius, int level, int numTotal, int b[])
	{	
		int mGetNum=total;
		int mNumTotal=numTotal;
		
		int i=level-1;
		int a[]=new int[level];
		
		for (int k = 0; k < a.length; k++) 
		{
			a[k] = k;
		}
		for(int m=0;m<b.length;m++)
		{
			a[m]=b[m];
		}
		
		
		if(mGetNum==mNumTotal)
		{
			
			return a;
		}
		
		
		while(a[i]==(i+level))
		{	
			i--;
		}
		a[i]=a[i]+1;
		for (int j = i + 1; j < level; j++) 
		{
			a[j] = a[i] + j - i;
		}
		return a;
	}
	
	public static double mean(double arr[][], int num)
	{
		double sum=0,mean=0;
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum=sum+arr[i][j];
			}
		}
		mean = sum/(num*num);
			
		return mean;
	}
	
	public static double variance(double arr[][], int num)
	{
		double sum=0,mean=0, var=0,sum_mean=0;
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum=sum+arr[i][j];
			}
		}
		mean = sum/(num*num);
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum_mean= sum_mean + ((arr[i][j]-mean)*(arr[i][j]-mean));
			}
		}
		return sum_mean;
	}
	
	public static double numerical_variance(double arr[][],double arr1[][],int num)
	{
		double sum=0,sum1=0,mean=0,mean1=0,sum_mean=0;
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum=sum+arr[i][j];
				sum1=sum1+arr1[i][j];
			}
		}
		mean = sum/(num*num);
		mean1=sum1/(num*num);
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum_mean= sum_mean + ((arr[i][j]-mean)*(arr1[i][j]-mean1));
			}
		}
		return sum_mean;
	}
	
	
}
