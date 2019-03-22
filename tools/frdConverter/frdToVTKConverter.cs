/*

Maciek Hawryłkiewicz (06.2015)
mamah(at)wp.pl


 */

using System;
using System.IO;
using System.Reflection;
using System.Collections.Generic;
using System.Text;

namespace mmh_ccx2vtk
{
	public partial class Program
	{

		// Komunikaty.cs -------------------------------------------------------------------------
		//INFORMATIVE TEXTS
		public static string TXT0000 = "*** CalculiX output file \".frd\" to vtk format converter by Maciej Marek Hawrylkiewicz mamah@wp.pl ***";

		public static string TXT0001 = ": Running on ";
		public static string TXT0002 = ": Read argument: ";
		public static string TXT0003 = ": Opening input file: \"{0}\"";
		public static string TXT0004 = ": Nodes read: {0}";
		public static string TXT0005 = ": Elements read: {0}";
		public static string TXT0006 = ": Reading result data from step: {0}";
		public static string TXT0007 = ": Nodal results block found:'{0}', defined number of result lines: {1}, defined number of items: {2}";
		public static string TXT0008 = ": Read {0} results lines in block: {1}";
		public static string TXT0009 = ": Read total number of steps: {0}";
		public static string TXT0010 = ": Writing new file: {0}";

		public static string TXT9999 = ": press Enter...";

		//ERROR MESSAGES
		public static string ERTXT0000 = "! (00) argument not found...";
		public static string ERTXT0001 = "! (01) Unexpected extension in file name \"{0}\", should be \".frd\"";
		public static string ERTXT0002 = "! (02) File not found...";

		public static string ERTXT0010 = "! (10) Cannot read the number of nodes in text \"{0}\",  in line: {1}";
		public static string ERTXT0011 = "! (11) Cannot read node data in line: {0}";
		public static string ERTXT0012 = "! (12) Declared number of nodes {0} is different than read node lines {1}";
		public static string ERTXT0013 = "! (13) Unexpected line code, in nodes block in line: {0}\"";
		public static string ERTXT0014 = "! (14) Cannot read element node number, in line: {0}";

		public static string ERTXT0020 = "! (20) Cannot read the number of elements in text \"{0}\",  in line: {1}";
		public static string ERTXT0021 = "! (21) Cannot read element node number, in line: {0}";
		public static string ERTXT0022 = "! (22) Unsupported element type in line: {0}";
		public static string ERTXT0023 = "! (23) Declared number of elements {0} is different than read element lines: {1}";
		public static string ERTXT0024 = "! (24) Unexpected line code in elements block, in line {0}";
		public static string ERTXT0025 = "! (25) The number of nodes numbers doesn't pass the element type, in line: {0}";

		public static string ERTXT0030 = "! (30) Cannot read the data block info in line {0}";
		public static string ERTXT0031 = "! (31) Can not read the dataset NAME or NUMBER of entities in line: {0}";
		public static string ERTXT0032 = "! (32) Can not read the dataset the entity NAME in line: {0}";
		public static string ERTXT0033 = "! (33) Can not read the results in line: {0}";
		public static string ERTXT0034 = "! (34) The number of read result lines {0} is different than declared number: {1}";

		public static void PrintHelp() {
			Console.WriteLine (":  ");
			Console.WriteLine (":  Usage:");
			Console.WriteLine (":  ");
			Console.WriteLine (":  mmh_ccx2vtk.exe [path]CalculiXoutputFileName.frd");
			Console.WriteLine (":  ");
			Console.WriteLine (":  The [path] string is optional, if not provided assembly will");
			Console.WriteLine (":  try to open specified \"frd\" file in current directory.");
			Console.WriteLine (":  ");
			Console.WriteLine (":  After succesful conversion the new file will be created using input file name ");
			Console.WriteLine (":  with suffix: \"_mmh.vtk\" in the input directory");
			Console.WriteLine (":  ");
			Console.WriteLine (":  Conversion supports all CalculiX element types and defined number of steps.");
			Console.WriteLine (":  The vector results are transalated as vectors, any other data are treat as scalar");
			Console.WriteLine (":  The name of Nodal Result Block in 'vtk' file and displayed in ParaView is as:");
			//Console.WriteLine (":      '[stepn_number]-block_name-item_name' where:");
			Console.WriteLine (":      'StepNumber_Blockname_ItemName' where:");
			Console.WriteLine (":  ");
//			Console.WriteLine (":            'step_number' is the user step number,");
//			Console.WriteLine (":            'block_name' is the name i.e. \"DISPR\"");
//			Console.WriteLine (":            'item_name' is the item name i.e. \"D1\"");
			Console.WriteLine (":            'StepNumber' is the user step number,");
			Console.WriteLine (":            'BlockName' is the name i.e. \"DISPR\"");
			Console.WriteLine (":            'ItemName' is the item name i.e. \"D1\"");
			Console.WriteLine (":  ");
			Console.WriteLine (":  It is the .NET/MONO assembly and it should work properly in Linux/Windows environment.");
			Console.WriteLine (":  ");
//			Console.WriteLine (":  Please inform me about bugs on mail: mamah@wp.pl");
			Console.WriteLine (":  ");
			Console.WriteLine (";) Have a nice step!");
			Console.WriteLine (":  ");
			Console.WriteLine (":  <<<MMH>>>");
			Console.WriteLine (":  ");
		}

		// Komunikaty.cs -------------------------------------------------------------------------


		// Struktury.cs -------------------------------------------------------------------------
		//node
		private struct tNode
		{
			public double X;
			public double Y;
			public double Z;
			public int vtk_nr;
			public bool IsFirst; //for element type 5-ccx 15 node wedge which seems to be not supported in VTK
		}

		//element
		private struct tElement
		{
			public int NR;
			public int TYP;
			public int[] N;
		}

		//the dictionary is really required here
		private struct tNodResBlRec {
			public double[] Data;
		}

		//Nodal Result Block
		private struct tNodResBlock {
			public int StepN; //step number
			public string Name; //name i.e. STRESSR
			public int Count; //number of lines with values
			public int DeclNumOfItems; //declared number of items
			public int NumOfItems; //as above but reduced when there are items supposed to be automatically calculated by cgx and there is no values in frd file
			public string[] ItemNameTab; //names of items (even with automated names i.e. 'ALL')
			public Dictionary<int, tNodResBlRec> DICT; //lines of values

			public tNodResBlock(string BlName, int TotNum, int ItNum, int StNum) {
				this.StepN = StNum;
				this.Name = BlName;
				this.Count = TotNum;
				this.DeclNumOfItems = ItNum;
				this.NumOfItems = ItNum;
				this.ItemNameTab = new string[ItNum];
				this.DICT = new Dictionary<int, tNodResBlRec> (TotNum);
			}
		}

		//------------ data containers ------------------------------------------------------
		private static Dictionary<int, tNode> NodeTab;
		private static tElement[] ElementTab;
		//result data blocks container
		private static List<tNodResBlock> DataTab;
		//------------------------------------------------------------------------------------

		// Struktury.cs -------------------------------------------------------------------------

		// Eksporty.cs -------------------------------------------------------------------------

		private static int VTKElemType(int Elind) {
			switch (ElementTab[Elind].TYP) {
			case 11:
				return 3;
			case 12:
				return 21;
			case 7:
				return 5;
			case 8:
				return 22;
			case 9:
				return 9;
			case 10:
				return 23;
			case 3:
				return 10;
			case 6:
				return 24;
			case 1:
				return 12;
			case 4:
				return 25;
			case 2:
			case 5:
				return 13;
			default:
				return 0;
			}
		}

		private static string VtkElementLine(int ElInd) {
			string tmpstr = "";
			int n, i;
			switch (ElementTab [ElInd].TYP) {
			case 11:
			case 12:
			case 7:
			case 8:
			case 9:
			case 10:
			case 3:
			case 6:
			case 1:

				n = ElementTab[ElInd].N.Length; tmpstr = string.Format("{0} ", n);
				for (i = 0; i < n; i++) 
					tmpstr += string.Format("{0, 8}", NodeTab[ElementTab[ElInd].N[i]].vtk_nr);
				break;
			case 4:			
				tmpstr = "20 ";
				for (i = 0; i < 12; i++) 
					tmpstr += string.Format("{0, 8}", NodeTab[ElementTab[ElInd].N[i]].vtk_nr);
				//last eight nodes has to be renumbered
				for (i = 16; i < 20; i++) 
					tmpstr += string.Format("{0, 8}", NodeTab[ElementTab[ElInd].N[i]].vtk_nr);
				for (i = 12; i < 16; i++) 
					tmpstr += string.Format("{0, 8}", NodeTab[ElementTab[ElInd].N[i]].vtk_nr);
				break;
			case 2:
			case 5:
				tmpstr = "6 ";
				tmpstr += string.Format("{0, 8}", NodeTab[ElementTab[ElInd].N[0]].vtk_nr);
				tmpstr += string.Format("{0, 8}", NodeTab[ElementTab[ElInd].N[2]].vtk_nr);
				tmpstr += string.Format("{0, 8}", NodeTab[ElementTab[ElInd].N[1]].vtk_nr);
				tmpstr += string.Format("{0, 8}", NodeTab[ElementTab[ElInd].N[3]].vtk_nr);
				tmpstr += string.Format("{0, 8}", NodeTab[ElementTab[ElInd].N[5]].vtk_nr);
				tmpstr += string.Format("{0, 8}", NodeTab[ElementTab[ElInd].N[4]].vtk_nr);
				break;
			}
			return tmpstr;
		}

		private static void MakeVTKFile() {

			//int TimeSteps = 10;
			int k = 0;

			Console.WriteLine ("Enter the number of Timesteps");
			int TimeSteps = Convert.ToInt32(Console.ReadLine());
			//TimeSteps = Console.ReadLine();
			//int TimeSteps = Convert.ToInt32(Console.ReadLine());

			for (k = 0; k < TimeSteps; k++){

			string vtk_plik = Path.GetDirectoryName(FrdFileName) + 
							  Path.DirectorySeparatorChar  + 
							  Path.GetFileNameWithoutExtension(FrdFileName) + "_mmh_" + k + ".vtk";

			string tmpstr;
			int i, j;

			Console.WriteLine (TXT0010, vtk_plik);

			// determine main nodes total number
			// elements ccx no.5 (15 node wedge) as they are not supported in VTK 
			// has to be exported as type ccx type 2 (VTK typ 13)
			// additional nodes has to be ommited
			int nn = 0;
	
			foreach (KeyValuePair<int, tNode> KVP in NodeTab)
				if (KVP.Value.IsFirst)
					nn++;
			// total number of elements
			int ne = ElementTab.Length;

			using (StreamWriter writer = new StreamWriter(vtk_plik)) {
				writer.WriteLine ("# vtk DataFile Version 3.0");
				writer.WriteLine ("MMH: CalculiX output '{0}' file to 'vtk' format", 
				                  Path.GetFileName(FrdFileName));
				writer.WriteLine ("ASCII");
				writer.WriteLine ("DATASET UNSTRUCTURED_GRID");
				writer.WriteLine ("POINTS {0} double", nn);

				//Nodes:
				List<int> tmpList = new List<int>(nn);
				foreach (KeyValuePair<int, tNode> KVP in NodeTab) 
					if (KVP.Value.IsFirst) {
						tmpList.Add (KVP.Key);
						tmpstr = string.Format (NFI, "{0, 16:0.00000000E+000}", KVP.Value.X) + 
								 string.Format (NFI, "{0, 16:0.00000000E+000}", KVP.Value.Y) +
								 string.Format (NFI, "{0, 16:0.00000000E+000}", KVP.Value.Z);
						writer.WriteLine (tmpstr);
					}

				// setting vtk numbers to nodes, 
				// now they will be numbered starting from 0
				// and counted as lines with coordinates
				tNode mN;
				for (i = 0; i < tmpList.Count; i++ ) {
					mN = NodeTab [tmpList [i]];
					mN.vtk_nr = i;
					NodeTab [tmpList [i]] = mN;
				}
				writer.WriteLine (); // tmpList.Clear ();

				// new number of nodes determine
				int totn = 0;
				for (i = 0; i < ElementTab.Length; i++)
					totn += (ElementTab[i].TYP == 5 ? 6 : ElementTab [i].N.Length);

				writer.WriteLine ("CELLS {0} {1}", ne, ne  + totn);

				for (i = 0; i < ne; i++) 
					writer.WriteLine (VtkElementLine(i));
				writer.WriteLine ();

				// element types:
				writer.WriteLine ("CELL_TYPES {0}", ne);
				for (i = 0; i < ne; i++)
					writer.WriteLine (VTKElemType(ElementTab[i].TYP));
				writer.WriteLine ();

				// cell data

				// not used in this case

				// point data
				writer.WriteLine ("POINT_DATA {0}", nn);

				//###################################################################
				//nodal results blocks iteration
				//for (i = 0; i < DataTab.Count; i++) {
				for (i = 0; i < DataTab.Count/TimeSteps; i++) {
					//writer.WriteLine ("DataTab.Count = {0}",  DataTab.Count);
				
				//	for (k = 0; k < 3; k++){

					if (DataTab [i].NumOfItems == 3) {
						//the blocks with 3 items are supposed to be a vectors
						//writer.WriteLine ("VECTORS [{0}]-{1} double", DataTab[i].StepN + 1,  DataTab[i].Name);
						writer.WriteLine ("VECTORS {1} double", DataTab[i + (k*(DataTab.Count/TimeSteps))].StepN + 1,  DataTab[i + (k*(DataTab.Count/TimeSteps))].Name);
						//writer.WriteLine ("DataTab[i].Name  = {0}",  DataTab[i].Name);
						//writer.WriteLine ("DataTab[i].NumOfItems  = {0}",  DataTab[i].NumOfItems);
						//line iteration
						foreach (KeyValuePair<int, tNodResBlRec > KVP in DataTab[i + (k*(DataTab.Count/TimeSteps))].DICT) {
							//check if the nodal data should be provided (element '5' case)
							if (NodeTab [KVP.Key].IsFirst) {
								tmpstr = "";
								//build line
								for (j = 0; j < DataTab[i + (k*(DataTab.Count/TimeSteps))].NumOfItems; j++) 
									tmpstr += string.Format(NFI, "{0, 16:0.0000000E+000}", KVP.Value.Data[j]);
								//write line
								writer.WriteLine (tmpstr);
							}
						}

					} else {
						// scalar results
						for (j = 0; j < DataTab[i + (k*(DataTab.Count/TimeSteps))].NumOfItems; j++) {
//							writer.WriteLine ("SCALARS [{0}]-{1}-{2} double", DataTab[i].StepN + 1, DataTab[i].Name, DataTab[i].ItemNameTab[j]);
							writer.WriteLine ("SCALARS {1}_{2} double", DataTab[i + (k*(DataTab.Count/TimeSteps))].StepN + 1, DataTab[i + (k*(DataTab.Count/TimeSteps))].Name, DataTab[i + (k*(DataTab.Count/TimeSteps))].ItemNameTab[j]);
							writer.WriteLine ("LOOKUP_TABLE default");
							//iteracja po liniach
							foreach (KeyValuePair<int, tNodResBlRec> KVP in DataTab[i + (k*(DataTab.Count/TimeSteps))].DICT) {
								if (NodeTab [KVP.Key].IsFirst) {
									tmpstr = string.Format (NFI, "{0, 16:0.0000000E+000}", KVP.Value.Data [j]);
									writer.WriteLine (tmpstr);
								}
							}
						}
					
					} //vector/scalar
				//	}
					// Write the VTK file must go here	- Kyle Davis 03/01/2019

				}  //block iteration
			}

			} //using
		} //method

		// Eksporty.cs -------------------------------------------------------------------------

		// Oblicyenia.cs -------------------------------------------------------------------------

				private static double FodX(double A, double B, double C, double x) {
			return x * x * x - A * x * x + B * x - C;
		}

		private static bool Poly3Root(double xa, double xb, double A, double B, double C, ref double x) {
			//x^3 - A * x^2 + B * x - C = 0
			bool result = false;

			double ya = FodX(A, B, C, xa);
			double yb = FodX(A, B, C, xb);

			double xc, yc;

			double D = 1e-6;

			if (Math.Abs (ya) < D) {
				result = true;
				x = xa;
			} else if (Math.Abs (yb) < D) {
				result = true;
				x = xb;
			} else if (ya > 0 && yb < 0 || ya < 0 && yb > 0) {
			
				for (;;) {

					xc = (xa + xb) / 2;
					yc = FodX (A, B, C, xc);

					if (Math.Abs (xb - xa) < D) {
						result = true;
						x = xc;
						break; 
					}

					if (ya < 0 && yc > 0 || ya > 0 && yc < 0) {
						xb = xc; yb = yc;
					} else {
						xa = xc; ya = yc;
					}

				}
			}
			return result;
		}

		// Oblicyenia.cs -------------------------------------------------------------------------


		//culture info
		private static System.Globalization.NumberFormatInfo NFI;

		//input file full name
		private static string FrdFileName = "";

		private static int StepNumber = -1;

		//nodes block
		private static int NumberOfNodes = 0;
		private static int NodeLineCounter = 0;

		//elements blok
		private static bool ELEM_FIRST_LINE = true;
		private static int NumberOfElems = 0;
		private static int CurElemType = 0;
		private static int ElemLineCounter = -1;	//starts from 0 to NumberOfElems - 1

		//nodal data block
		private static int DataBkInd = -1;
		private static int DataLineCounter = 0;

		private static bool NameOfInputFile(ref string mPlik)
		{
			string mpath = "";
			string mext = "";
			//check if file name has extension
			mext = Path.GetExtension(mPlik);
			if (string.IsNullOrEmpty(mext)) {
				//add proper extension if it was not given by user
				mPlik += ".frd";
			} else if (mext != ".frd") {
				//given improper extension 
				Console.WriteLine(ERTXT0001, mext);
				Console.Write(TXT9999);
				Console.ReadLine();
				return false;
			}
			//read path if it was given:
			mpath = Path.GetDirectoryName(mPlik);
			//if there is no path in argument then use assembly execution path:
			if (string.IsNullOrEmpty(mpath))
				mpath = Path.GetDirectoryName(Assembly.GetExecutingAssembly().CodeBase);
			//full path creation
			mPlik = mpath + Path.DirectorySeparatorChar + Path.GetFileNameWithoutExtension(mPlik) + mext;
			//removing prefix "file:" Unix / "file:\" Windows
			if (mPlik.Contains("file:")) {
				if (System.Environment.OSVersion.VersionString.Contains("Unix")) {				
					mPlik = mPlik.Substring(5);
				} else {
					mPlik = mPlik.Substring(6);
				}
			}

			Console.WriteLine(TXT0003, mPlik);

			return (File.Exists(mPlik));
		}

		private static bool ProcessNodeHeadLine(ref string mLinia, ref int NrL) {
			string[] Tstr = mLinia.Split(new Char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
			//number of nodes:
			try {
				NumberOfNodes = Convert.ToInt32(Tstr[1]);
			} catch {
				Console.WriteLine(ERTXT0010, Tstr[1], NrL);
				return false;
			}
			//set node database
			NodeTab = new System.Collections.Generic.Dictionary<int, tNode>(NumberOfNodes);
			return true;
		}

		private static bool ProcessNodeLine(ref string mLinia, ref int NrL, ref bool KONIEC)
		{
			tNode mNode = default(tNode);
			int tmpNR = 0;

			if (mLinia.IndexOf("-1") == 1) {
				try {
					NodeLineCounter += 1;
					tmpNR = Convert.ToInt32(mLinia.Substring(3, 10), NFI);
					mNode.X = Convert.ToDouble(mLinia.Substring(13, 12), NFI);
					mNode.Y = Convert.ToDouble(mLinia.Substring(25, 12), NFI);
					mNode.Z = Convert.ToDouble(mLinia.Substring(37, 12), NFI);
					mNode.IsFirst = true;
					NodeTab.Add(tmpNR, mNode);
				} catch {
					Console.WriteLine(ERTXT0011, NrL);
					return false;
				}
			} else if (mLinia.IndexOf("-3") == 1) {
				//end of reading
				if (NodeTab.Count == NodeLineCounter) {
					//read number is correct
					KONIEC = true;
					Console.WriteLine(TXT0004, NodeLineCounter);
				} else {
					//the numbers differs
					Console.WriteLine(ERTXT0012, NodeLineCounter, NodeTab.Count);
					return false;
				}
			} else {
				Console.WriteLine(ERTXT0013, NrL);
				return false;
			}

			return true;
		}

		private static int NumOfElemNodes(int ElementType) {

			switch (ElementType) {
			case 11:
				return 2;
			case 12:
			case 7:
				return 3;
			case 8:
				return 6;
			case 9:
				return 4;
			case 10:
				return 8;
			case 3:
				return 4;
			case 6:
				return 10;
			case 1:
				return 8;
			case 4:
				return 20;
			case 2:
				return 6;
			case 5:
				return 15;
			default:
				return 0;
			}

		}

		private static bool ReadNodesListLine(ref string[] tmpstr, int NrL, int ExpNum) {
			int tmpNr;
			if (tmpstr.Length == ExpNum + 1) {
				for (int i = 0; i < ExpNum; i++) {
					try {
						tmpNr = Convert.ToInt32 (tmpstr [i + 1]);
						ElementTab [ElemLineCounter].N [i] = tmpNr;
					} catch {
						Console.WriteLine(ERTXT0014, NrL.ToString());
						return false;
					}
				}
			} else {
				Console.WriteLine(ERTXT0025, NrL.ToString());
				return false;
			} 
			return true;
		}

		private static bool ProcessElementHeadLine(ref string MLinia, ref int NrL) {
			string[] TStr = MLinia.Split(new Char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
			try {
				NumberOfElems = Convert.ToInt32(TStr[1]);
			} catch {
				Console.WriteLine(ERTXT0020, TStr[1], NrL.ToString());
				return false;
			}
			ElementTab = new tElement[NumberOfElems] ;
			return true;
		}

		private static bool ProcessElementLine(ref string MLinia, ref int NrL, ref bool KONIEC)	{

			string[] TStr = null;
			int tmpNR = 0;
			int i;

			if (MLinia.IndexOf("-1") == 1) {

				ElemLineCounter++;
				TStr = MLinia.Split(new Char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);

				ElementTab[ElemLineCounter].NR = Convert.ToInt32(TStr[1]);
				CurElemType = Convert.ToInt32(TStr[2]);
				ElementTab[ElemLineCounter].TYP = CurElemType;

				ElementTab [ElemLineCounter].N = new int[ NumOfElemNodes (CurElemType) ];

			} else if (MLinia.IndexOf("-2") == 1) {
				TStr = MLinia.Split(new Char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);

				switch (CurElemType) {
				case 11:
					if (!ReadNodesListLine(ref TStr,NrL,2)) return false;
					break;
				case 12:
				case 7:
					if (!ReadNodesListLine(ref TStr,NrL,3)) return false;
					break;
				case 2:
				case 8:
					if (!ReadNodesListLine(ref TStr,NrL,6)) return false;
					break;
				case 3:
				case 9:
					if (!ReadNodesListLine(ref TStr,NrL,4)) return false;
					break;
				case 1:
				case 10:
					if (!ReadNodesListLine(ref TStr,NrL,8)) return false;
					break;
				case 6:
					if (!ReadNodesListLine(ref TStr,NrL,10)) return false;
					break;
				case 4:
					// two lines with ten nodes 
					if (TStr.Length == 11) {
						int offset = (ELEM_FIRST_LINE ? 0 : 10);
						for (i = 0; i < 10; i++) {
							try {
								tmpNR = Convert.ToInt32 (TStr [i + 1]);
								ElementTab [ElemLineCounter].N [i + offset] = tmpNR;

							} catch {
								Console.WriteLine (ERTXT0021, NrL.ToString ()); 
								return false;
							}
						}
					} else {
						Console.WriteLine (ERTXT0025, NrL.ToString ());
						return false;
					}
					ELEM_FIRST_LINE = !ELEM_FIRST_LINE;
					break;
				case 5:
					// two line 10 + 5 nodes - element not supported in VTK - will be processed as element type 2 (6 in VTK)
					if (ELEM_FIRST_LINE) {
						if (TStr.Length == 11) {

							for ( i = 0; i < 10; i++ ) {
								try {
									tmpNR = Convert.ToInt32(TStr[i + 1]);
									ElementTab[ElemLineCounter].N[i] = tmpNR;
								} catch {
									Console.WriteLine(ERTXT0021, NrL.ToString());
									return false;
								}
							}
						} else {
							Console.WriteLine(ERTXT0025, NrL.ToString());
							return false;
						}
						ELEM_FIRST_LINE = false;
					} else {
						//last five nodes
						if (TStr.Length == 6)
							for ( i = 0; i < 5; i++ ) 
								try {
									tmpNR = Convert.ToInt32(TStr[i + 1]);
									ElementTab[ElemLineCounter].N[i + 10] = tmpNR;
								} catch {
									Console.WriteLine(ERTXT0021, NrL.ToString());
									return false;
								}
						else {
							Console.WriteLine(ERTXT0025, NrL.ToString());
							return false;
						}
						// only first six nodes 0..5 will be transfered to VTK
						tNode mNode;
						for (i = 6; i < 15; i++ ) {
							mNode = NodeTab[ElementTab[ElemLineCounter].N[i]];
							mNode.IsFirst = false;
							NodeTab [ElementTab [ElemLineCounter].N [i]] = mNode;
						}

						ELEM_FIRST_LINE = true;
					}
					break;
				default: {
					Console.WriteLine(ERTXT0022, NrL.ToString());
					return false;
					}

				}

			} else if (MLinia.IndexOf("-3") == 1) {
				//end of elements block
				if (ElemLineCounter == NumberOfElems - 1) {
					KONIEC = true;
					Console.WriteLine(TXT0005, NumberOfElems.ToString());
				} else {
					Console.WriteLine(ERTXT0023, NumberOfElems, ElemLineCounter);
					return false;
				}
			} else {
				//improper line starting code
				Console.WriteLine(ERTXT0024, NrL.ToString());
				return false;
			}

			return true;
		}
	
		private static bool ProcessDataLine(ref string MLinia, ref int NrL, ref bool KONIEC) {
			tNodResBlRec mNodResBlRec = new tNodResBlRec();
			int tmpNr, k = 13;
			if (MLinia.IndexOf ("-1") == 1) {
				DataLineCounter++;
				mNodResBlRec.Data = new double[DataTab[DataBkInd].NumOfItems];
				try {
					tmpNr = Convert.ToInt32(MLinia.Substring(3, 10));
					for ( int i = 0; i < DataTab[DataBkInd].NumOfItems; i++ ) {
						mNodResBlRec.Data[i] = Convert.ToDouble(MLinia.Substring(k, 12), NFI);
						k += 12;
					}
					DataTab[DataBkInd].DICT.Add(tmpNr, mNodResBlRec);
				} catch {
					Console.WriteLine(ERTXT0033, NrL.ToString());
					return false;
				}
			} else if (MLinia.IndexOf("-3") == 1) {
				if (DataTab [DataBkInd].Count == DataLineCounter) {
					KONIEC = true;
					Console.WriteLine (TXT0008, DataLineCounter.ToString(), DataTab[DataBkInd].Name);
					DataLineCounter = 0;
				} else {
					Console.WriteLine (ERTXT0034, DataLineCounter, DataTab[DataBkInd].Count);
					return false;
				}
			}
			return true;
		}

		private static bool ReadFrdData()
		{
			bool CZYT_WEZLY = false;
			bool CZYT_ELEMENTY = false;
			bool CZYT_DANE = false;

			bool OK_WEZLY = false;
			bool OK_ELEMENTY = false;

			bool OK_BLOK = false;

			int NrL = 0;
			string mLinia = "";

			//linia "100C"
			int NumOfNodeData = 0;
			//linia "-4"
			string DataName = "";
			int NumOfEnt = 0;

			//linia "-5"
			int EntNamesInd = -1;

			tNodResBlock mNodResBlock;

			DataTab = new List<tNodResBlock> ();

			using (StreamReader reader = File.OpenText(FrdFileName)) {

				while (reader.Peek() >= 0) {
					mLinia = reader.ReadLine();	NrL++;

					//FLAGS
					if (mLinia.IndexOf ("2C") == 4) {
						if ( ProcessNodeHeadLine(ref mLinia, ref NrL) ) {
							CZYT_WEZLY = true;
							continue;
						} else 
							break;
					}

					if (mLinia.IndexOf ("3C") == 4) {
						if ( ProcessElementHeadLine (ref mLinia, ref NrL) ) {
							CZYT_ELEMENTY = true;
							continue;
						} else 
							break;
					}

//					//Parameter Header Record:
//					if (mLinia.IndexOf ("    1P") == 0)
//						continue;

					//Nodal Result Block:
					if (mLinia.IndexOf("100C") == 2) {
						/*
						Format:(1X,’ 100’,’C’,6A1,E12.5,I12,20A1,I2,I5,10A1,I2)
						Values: KEY,CODE,SETNAME,VALUE,NUMNOD,TEXT,ICTYPE,NUMSTP,ANALYS,FORMAT
						 */
						try {

							NumOfNodeData = Convert.ToInt32(mLinia.Substring(24, 12));
							int SNr = Convert.ToInt32 (mLinia.Substring(58, 5)) - 1;

							if (StepNumber != SNr) {
								StepNumber = SNr;
								Console.WriteLine(TXT0006, StepNumber + 1);
							}
						} catch {
							Console.WriteLine (ERTXT0030, NrL);
							return false;
						}
						continue;
					}

					/*
							Format:(1X,I2,2X,8A1,2I5)
							Values: KEY, NAME, NCOMPS, IRTYPE
							Where: KEY = -4
							NAME = Dataset name to be used in the menu
							NCOMPS = Number of entities
							IRTYPE = 1 Nodal data, material independent
							2 Nodal data, material dependant
							3 Element data at nodes (not used)
					 */
					if (mLinia.IndexOf ("-4") == 1) {
						try {
//							DataName = mLinia.Substring(5, 8).Trim();
							DataName = mLinia.Substring(5, 8).Replace(" ", "");
							NumOfEnt = Convert.ToInt32 (mLinia.Substring(13, 5));
						} catch {
							Console.WriteLine (ERTXT0031, NrL);
							return false;
						}


						Console.WriteLine(TXT0007, DataName, NumOfNodeData, NumOfEnt);

						mNodResBlock = new tNodResBlock (DataName, NumOfNodeData, NumOfEnt, StepNumber);
						DataTab.Add(mNodResBlock);
						DataBkInd++; //current DataTab index (current block of data index)
						EntNamesInd = -1; //prepare to read items from lines code '-5'
						continue;
					}

					if (mLinia.IndexOf ("-5") == 1) {
						EntNamesInd++;
						try {
//							DataTab[DataBkInd].ItemNameTab[EntNamesInd] = mLinia.Substring(5,8).Trim();
							DataTab[DataBkInd].ItemNameTab[EntNamesInd] = mLinia.Substring(5,8).Replace(" ", "");
						} catch {
							Console.WriteLine (ERTXT0032, NrL);
							return false;
						}

						// if there is an item to be auto calculated by cgx 
						// reduce number of expected values in data lines
						try {
							if ( Convert.ToInt32 (mLinia.Substring(33, 5)) == 1 ) {
								mNodResBlock = DataTab[DataBkInd];
								mNodResBlock.NumOfItems--;
								DataTab[DataBkInd] = mNodResBlock;
							}
						} catch {/*nothing to do*/}	

						//check if this is the last code '-5' line
						if (DataTab [DataBkInd].DeclNumOfItems == EntNamesInd + 1)
							CZYT_DANE = true;

						continue;
					}

					//Processing
					if (CZYT_WEZLY) {
						if (ProcessNodeLine (ref mLinia, ref NrL, ref OK_WEZLY)) {
							if (OK_WEZLY)
								CZYT_WEZLY = false;
						} else 
							break;

					} else if (OK_WEZLY && CZYT_ELEMENTY) {
						if (ProcessElementLine (ref mLinia, ref NrL, ref OK_ELEMENTY)) {
							if (OK_ELEMENTY)
								CZYT_ELEMENTY = false;
						} else 
							break;

					} else if (OK_WEZLY && OK_ELEMENTY && CZYT_DANE) {
						if (ProcessDataLine (ref mLinia, ref NrL, ref OK_BLOK)) {
							if (OK_BLOK) {
								OK_BLOK = false;
								CZYT_DANE = false;
								DataLineCounter = 0;
							}
						} else
							break;
					} //else
				} //while
			} //using
			//	nodal result datas are not required, the model is enough
			return (OK_WEZLY & OK_ELEMENTY);
		}

		public static void Main(string[] args)
		{

			NFI = new System.Globalization.NumberFormatInfo ();
			NFI.NumberDecimalSeparator = ".";

			Console.WriteLine (TXT0000);
			Console.WriteLine (TXT0001 + System.Environment.OSVersion.VersionString);

			//for debugging
			//args = new string[] { "/media/mmh/E_Windows/_CALCULIX_BCONV_WORKING/dupa/dupa.frd" };

			if (args.Length != 1) {
				Console.WriteLine(ERTXT0000);
				PrintHelp ();
			} else {

				Console.WriteLine (TXT0002 + args[0]);
				FrdFileName = args[0];

				if (NameOfInputFile(ref FrdFileName)) {
					if (ReadFrdData()) {
						Console.WriteLine (TXT0009, StepNumber + 1);
						MakeVTKFile ();
					}
				} else {
					Console.WriteLine(ERTXT0002);
					PrintHelp ();
				}
			}

			Console.WriteLine();
			//Console.Write(TXT9999);
			//Console.ReadLine();

		}
	}
}
