
/*

Repair -- an implementation of Larsson and Moffat's compression and
decompression algorithms.
Copyright (C) 2010-current_year Gonzalo Navarro

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Author's contact: Gonzalo Navarro, Dept. of Computer Science, University of
Chile. Blanco Encalada 2120, Santiago, Chile. gnavarro@dcc.uchile.cl

*/

int PRNC = 0;  // print current sequence C (verbose!)
int PRNR = 0;  // print active pairs in the heap (verbose!)
int PRNP = 0;  // print forming pairs
int PRNL = 0;  // print progress on text scan

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "basics.h"
#include "records.h"
#include "hash.h"
#include "heap.h"

float factor = 0.75; // 1/extra space overhead; set closer to 1 for smaller and
		     // slower execution
int minsize = 256; // to avoid many reallocs at small sizes, should be ok as is

relong u; // |text| and later current |C| with gaps

unsigned short *sC; // shorts version of C

relong *C; // compressed text, pero tiene -ptrs al mismo texto

relong c;  // real |C|

int alph; // max used terminal symbol

int n; // |R|

Tlist *L; // |L| = c;

Thash Hash; // hash table of pairs

Theap Heap; // special heap of pairs

Trarray Rec; // records

int chars[256];

char map[256];

unsigned char *text; // first pass over this char array

int MB;

int did0 = 0; // did prepare0 or not

	// first pass, inplace on text, up to completing 256 symbols

void prepare0 (relong len)

  { relong i;
    int id;
    Tpair pair;
    alph = 0;
    for (i=0;i<256;i++) chars[i] = -1;
    for (i=0;i<len;i++) 
	{ unsigned char x = text[i];
	  if (chars[x] == -1) chars[x] = alph++;
	  text[i] = chars[x];
	}
    for (i=0;i<256;i++) 
	if (chars[i] != -1) map[chars[i]] = i;
    n = alph;
    Rec = createRecords(factor,minsize);
    Heap = createHeap(len,&Rec,factor,minsize);
    Hash = createHash(256*256,&Rec);
    assocRecords (&Rec,&Hash,&Heap,NULL);
    if ((len /1024/1024) * 3 * sizeof(relong) <= MB) return;
    did0 = 1;
    for (i=0;i<len-1;i++) 
	{ pair.left = text[i]; pair.right = text[i+1];
	  id = searchHash (Hash,pair);
    	  if (id == -1) // new pair, insert
	     { id = insertRecord (&Rec,pair);
	     }
	  else 
	     { incFreq (&Heap,id);
	     }
if (PRNL && (i%1000000 == 0)) printf ("Processed %lli chars\n",i);
	}
    purgeHeap (&Heap);
  }

	// second pass, convert text to shorts

void prepare1 (relong len)

  { relong i;
    if ((len /1024/1024) * 3 * sizeof(relong) <= MB)
       { sC = NULL; return; }
    sC = malloc (len * sizeof(short));
    for (i = 0; i < len; i++) sC[i] = text[i];
    free (text);
  }

	// third pass, use the full machinery

void prepare (relong len)

  { relong i;
    int id;
    Tpair pair;
    c = u = len;
    C = (void*)malloc(u*sizeof(relong));
    if (sC != NULL)
       { for (i=0;i<u;i++) C[i] = sC[i];
         free (sC);
       }
    else
       { for (i=0;i<u;i++) C[i] = text[i];
         free (text);
       }
    L = (void*)malloc(u*sizeof(Tlist));
    assocRecords (&Rec,&Hash,&Heap,L);
    for (i=0;i<c-1;i++) 
	{ pair.left = C[i]; pair.right = C[i+1];
	  id = searchHash (Hash,pair);
    	  if (id == -1) // non existing => occurs only once, don't create
	     { if (did0) L[i].prev = NullFreq;
	       else { id = insertRecord (&Rec,pair);
	              L[i].prev = -id-1; 
	              Rec.records[id].cpos = i;
	            } 
	       L[i].next = -1;
	     }
	  else 
	     { if (Rec.records[id].cpos == -1) // first time I see it this pass
	          { L[i].next = -1;
	          }
	       else
	          { L[i].next = Rec.records[id].cpos;
	            L[L[i].next].prev = i;
	            if (!did0) incFreq (&Heap,id);
	          }
	       L[i].prev = -id-1; 
	       Rec.records[id].cpos = i;
	     }
if (PRNL && (i%1000000 == 0)) printf ("Processed %lli chars\n",i);
	}
    L[i].prev = NullFreq; L[i].next = -1;
    purgeHeap (&Heap);
  }

void prnSym(int c)
   { if (c < alph) printf("%c",map[c]); else printf ("%i",c);
   }

void prntext (relong len)

  { relong i = 0;
    printf ("C[1..%lli] = ",len);
    while (i<len)
      { prnSym(text[i]);
  	printf (" ");
	i++; 
      }
    printf ("\n\n");
  }

void prnsC (relong len)

  { relong i = 0;
    printf ("C[1..%lli] = ",len);
    while (i<len)
      { prnSym(sC[i]);
  	printf (" ");
	i++; 
      }
    printf ("\n\n");
  }

void prnC (void)

  { relong i = 0;
    printf ("C[1..%lli] = ",c);
    while (i<u)
      { prnSym((int)C[i]);
  	printf (" ");
	i++; 
	if ((i<u) && (C[i] < 0)) i = -C[i]-1; 
      }
    printf ("\n\n");
  }

void prnRec (void)

  { int i;
    printf ("Active pairs:\n");
    for (i=0;i<Rec.size;i++)
        { printf ("\t(");
	  prnSym(Rec.records[i].pair.left);
          printf (",");
	  prnSym(Rec.records[i].pair.right);
	  printf ("), %lli occs\n", Rec.records[i].freq);
	}
    printf ("\n");
  }

	// first pass, work directly on the text chars up to 256 symbols

relong repair0 (relong len, FILE *R)

  { int oid,id;
    relong cpos,pos;
    Trecord *rec,*orec;
    Tpair pair;
    int left,right;
    if (fwrite(&alph,sizeof(int),1,R) != 1) return -1;
    if (fwrite(map,sizeof(char),alph,R) != alph) return -1;
if (PRNL) printf ("--- first stage, n = %lli\n",len);
if (PRNC) prntext(len); 
    while (n < 256)
      { 
        if ((len /1024/1024) * 3 * sizeof(relong) <= MB) return len;
	else if (PRNP) printf ("Avoiding to use %lli MB\n",
			       (len /1024/1024) * 3 * sizeof(relong));
if (PRNR) prnRec();
	oid = extractMax(&Heap);
	if (oid == -1) break; // the end!!
 	orec = &Rec.records[oid];
	if (fwrite (&orec->pair,sizeof(Tpair),1,R) != 1) return -1;
        left = orec->pair.left; right = orec->pair.right;
if (PRNP) 
    { printf("Chosen pair %i = (",n);
      prnSym(orec->pair.left);
      printf(",");
      prnSym(orec->pair.right);
      printf(") (%lli occs)\n",orec->freq);
    }
	pos = 0;
	for (cpos=0;cpos<len-1;cpos++)
	   { if ((text[cpos] != left) || (text[cpos+1] != right))
		text[pos] = text[cpos];
	     else // occurrence of the pair to modify
		{   // decr freqs of pairs that disappear
		  if (pos > 0) 
		     { pair.left = text[pos-1]; pair.right = text[cpos];
		       id = searchHash(Hash,pair);
		       if (id != -1) // may not exist if purgeHeap'd
	                  { if (id != oid) decFreq(&Heap,id); //not to my pair!
			  }
		     }
		  if (cpos < len-2) 
		     { pair.left = text[cpos+1]; pair.right = text[cpos+2];
		       id = searchHash(Hash,pair);
		       if (id != -1) // may not exist if purgeHeap'd
	                  { if (id != oid) decFreq(&Heap,id); //not to my pair!
			  }
		     }
			// create or incr freqs of pairs that appear
		  if (pos > 0) 
		     { pair.left = text[pos-1];
		       pair.right = n;
		       id = searchHash(Hash,pair);
	               if (id == -1) // new pair, insert
		          { id = insertRecord (&Rec,pair);
		          }
	               else 
 		          { incFreq (&Heap,id);
	                  }
	             }
		  if (cpos < len-2) 
		     { pair.left = n; pair.right = text[cpos+2];
		       id = searchHash(Hash,pair);
	               if (id == -1) // new pair, insert
		          { id = insertRecord (&Rec,pair);
		          }
	               else 
 		          { incFreq (&Heap,id);
	                  }
	             }
			// actually do the substitution
		  text[pos] = n; cpos++;
		}
	      pos++;
	   }
         if (cpos == len-1) text[pos++] = text[cpos];
	 len = pos;
if (PRNC) prntext(len);
	 removeRecord (&Rec,oid);
	 n++;
	 purgeHeap(&Heap); // remove freq 1 from heap
         text = (void*)realloc ((void*)text,len*sizeof(char)); // should not move it
       }
     purgeHeap(&Heap); // remove freq 1 from heap, if it exited for oid=-1
     return len;
   }

	// second pass, work on a version of shorts, until the structure fits
	// in MB, or until shorts are insufficient

relong repair1 (relong len, FILE *R)

  { int oid,id;
    relong cpos,pos,ofreq;
    Trecord *rec,*orec;
    Tpair pair;
    int left,right;
if (PRNL) printf ("--- second stage, n=%lli\n",len);
if (PRNC) prnsC(len);
    while (n < 1 << 16)
      { 
        if ((len /1024/1024) * 3 * sizeof(relong) <= MB) return len;
	else if (PRNP) printf ("Avoiding to use %lli MB\n",
			       (len /1024/1024) * 3 * sizeof(relong));
if (PRNR) prnRec();
	oid = extractMax(&Heap);
	if (oid == -1) break; // the end!!
 	orec = &Rec.records[oid];
	if (fwrite (&orec->pair,sizeof(Tpair),1,R) != 1) return -1;
        left = orec->pair.left; right = orec->pair.right;
        ofreq = orec->freq;
if (PRNP) 
    { printf("Chosen pair %i = (",n);
      prnSym(orec->pair.left);
      printf(",");
      prnSym(orec->pair.right);
      printf(") (%lli occs)\n",orec->freq);
    }
	pos = 0;
	for (cpos=0;cpos<len-1;cpos++)
	   { if ((sC[cpos] != left) || (sC[cpos+1] != right))
		sC[pos] = sC[cpos];
	     else // occurrence of the pair to modify
		{   // decr freqs of pairs that disappear
		  if (pos > 0) 
		     { pair.left = sC[pos-1]; pair.right = sC[cpos];
		       id = searchHash(Hash,pair);
		       if (id != -1) // may not exist if purgeHeap'd
	                  { if (id != oid) decFreq(&Heap,id); //not to my pair!
			  }
		     }
		  if (cpos < len-2) 
		     { pair.left = sC[cpos+1]; pair.right = sC[cpos+2];
		       id = searchHash(Hash,pair);
		       if (id != -1) // may not exist if purgeHeap'd
	                  { if (id != oid) decFreq(&Heap,id); //not to my pair!
			  }
		     }
			// create or incr freqs of pairs that appear
		  if (pos > 0) 
		     { pair.left = sC[pos-1];
		       pair.right = n;
		       id = searchHash(Hash,pair);
	               if (id == -1) // new pair, insert
		          { id = insertRecord (&Rec,pair);
		          }
	               else 
 		          { incFreq (&Heap,id);
	                  }
	             }
		  if (cpos < len-2) 
		     { pair.left = n; pair.right = sC[cpos+2];
		       id = searchHash(Hash,pair);
	               if (id == -1) // new pair, insert
		          { id = insertRecord (&Rec,pair);
		          }
	               else 
 		          { incFreq (&Heap,id);
	                  }
	             }
			// actually do the substitution
		  sC[pos] = n; cpos++;
		}
	      pos++;
	   }
         if (cpos == len-1) sC[pos++] = sC[cpos];
	 len = pos;
if (PRNC) prnsC(len);
	 removeRecord (&Rec,oid);
	 n++;
	 purgeHeap(&Heap); // remove freq 1 from heap
         sC = realloc (sC,len*sizeof(short)); // should not move it
       }
     purgeHeap(&Heap); // remove freq 1 from heap, if it exited for oid=-1
     return len;
   }

	// to avoid another list to access the sparse C we thread it using the
	// empty space. if next cell of an active cell is negative, it is 
	// (minus) a ptr to the next occ. idem previous cell to previous occ,
	// except that next ptr dominates over prev ptr if they must be in
 	// the same cell. but in this case one can find prev in O(1) anyway.

relong repair (FILE *R)

  { int oid,id;
    relong cpos;
    Trecord *rec,*orec;
    Tpair pair;
if (PRNL) printf ("--- third stage, n=%lli\n",c);
if (PRNC) prnC();
    while (n+1 > 0)
      { 
if (PRNR) prnRec();
	oid = extractMax(&Heap);
	if (oid == -1) break; // the end!!
 	orec = &Rec.records[oid];
 	cpos = orec->cpos;
	if (fwrite (&orec->pair,sizeof(Tpair),1,R) != 1) return -1;
if (PRNP) 
    { printf("Chosen pair %i = (",n);
      prnSym(orec->pair.left);
      printf(",");
      prnSym(orec->pair.right);
      printf(") (%lli occs)\n",orec->freq);
    }
	while (cpos != -1)
	   { relong ant,sgte,ssgte; 
		// replacing bc->e in abcd, b = cpos, c = sgte, d = ssgte
	     if (C[cpos+1] < 0) sgte = -C[cpos+1]-1; 
	     else sgte = cpos+1; 
	     if ((sgte+1 < u) && (C[sgte+1] < 0)) ssgte = -C[sgte+1]-1;
	     else ssgte = sgte+1; 
		// remove bc from L
	     if (L[cpos].next != -1) L[L[cpos].next].prev = -oid-1;
	     orec->cpos = L[cpos].next;
	     if (ssgte != u) // there is ssgte
		{ 	// remove occ of cd
		  pair.left = C[sgte]; pair.right = C[ssgte];
		  id = searchHash(Hash,pair);
		  if (id != -1) // may not exist if purgeHeap'd
	             { if (id != oid) decFreq (&Heap,id); // not to my pair!
		       if (L[sgte].prev != NullFreq) //still exists(not removed)
 		          { rec = &Rec.records[id];
		            if (L[sgte].prev < 0) // this cd is head of its list
		               rec->cpos = L[sgte].next;
		            else L[L[sgte].prev].next = L[sgte].next;
		            if (L[sgte].next != -1) // not tail of its list
		               L[L[sgte].next].prev = L[sgte].prev;
			  }
		     }
			// create occ of ed
	          pair.left = n;
		  id = searchHash(Hash,pair);
	          if (id == -1) // new pair, insert
		     { id = insertRecord (&Rec,pair);
 		       rec = &Rec.records[id];
	               L[cpos].next = -1;
		     }
	          else 
 		     { incFreq (&Heap,id);
 		       rec = &Rec.records[id]; 
	               L[cpos].next = rec->cpos;
	               L[L[cpos].next].prev = cpos;
	             }
	          L[cpos].prev = -id-1;
	          rec->cpos = cpos;
		}
	     if (cpos != 0) // there is ant
		{ 	// remove occ of ab
	          if (C[cpos-1] < 0) 
		     { ant = -C[cpos-1]-1; 
		       if (ant == cpos) // sgte and ant clashed -> 1 hole
			  ant = cpos-2;
		     }
	          else ant = cpos-1; 
		  pair.left = C[ant]; pair.right = C[cpos];
		  id = searchHash(Hash,pair);
		  if (id != -1) // may not exist if purgeHeap'd
	             { if (id != oid) decFreq (&Heap,id); // not to my pair!
		       if (L[ant].prev != NullFreq) //still exists (not removed)
 		          { rec = &Rec.records[id];
		            if (L[ant].prev < 0) // this ab is head of its list
		                 rec->cpos = L[ant].next;
		            else L[L[ant].prev].next = L[ant].next;
		            if (L[ant].next != -1) // it is not tail of its list
		               L[L[ant].next].prev = L[ant].prev;
			  }
		     }
			// create occ of ae
	          pair.right = n;
		  id = searchHash(Hash,pair);
	          if (id == -1) // new pair, insert
		     { id = insertRecord(&Rec,pair);
 		       rec = &Rec.records[id];
	               L[ant].next = -1;
	             }
	          else 
	             { incFreq (&Heap,id);
 		       rec = &Rec.records[id];
	               L[ant].next = rec->cpos;
	               L[L[ant].next].prev = ant;
	             }
	          L[ant].prev = -id-1;
	          rec->cpos = ant;
		}
	     C[cpos] = n;
	     if (ssgte != u) C[ssgte-1] = -cpos-1;
	     C[cpos+1] = -ssgte-1;
	     c--;
 	     orec = &Rec.records[oid]; // just in case of Rec.records realloc'd
	     cpos = orec->cpos;
	   }
if (PRNC) prnC();
	 removeRecord (&Rec,oid);
	 n++;
	 purgeHeap(&Heap); // remove freq 1 from heap
	 if (c < factor * u) // compact C
	    { relong i,ni;
	      i = 0;
	      for (ni=0;ni<c-1;ni++) 
		{ C[ni] = C[i];
		  L[ni] = L[i];
		  if (L[ni].prev < 0) 
		     { if (L[ni].prev != NullFreq) // real ptr
			  Rec.records[-L[ni].prev-1].cpos = ni; 
		     }
		  else L[L[ni].prev].next = ni;
		  if (L[ni].next != -1) L[L[ni].next].prev = ni;
		  i++; if (C[i] < 0) i = -C[i]-1;
		}
	      C[ni] = C[i];
	      u = c;
	      C = realloc (C, c * sizeof(relong));
	      L = realloc (L, c * sizeof(Tlist));
              assocRecords (&Rec,&Hash,&Heap,L);
	    }
       }
     return 0;
   }

void main (int argc, char **argv)

   { char fname[1024];
     FILE *Tf,*Rf,*Cf;
     relong i,olen,len;
     struct stat s;
     if (argc != 3)
	{ fprintf (stderr,"Usage: %s <filename> <MB>\n"
			  "Compresses <filename> with repair and creates "
			  "<filename>.ext compressed files\n"
			  "Use <MB> to give an approximate idea of how much "
			  "main memory you have for compression\n",argv[0]);
	  exit(1);
	}
     MB = atoi(argv[2]);
     if (MB <= 0)
	{ fprintf (stderr,"Usage: %s <filename> <MB>\n"
			  "Compresses <filename> with repair and creates "
			  "<filename>.ext compressed files\n"
			  "Use <MB> to give an approximate idea of how much "
			  "main memory you have for compression\n",argv[0]);
	  exit(1);
	}
     if (stat (argv[1],&s) != 0)
	{ fprintf (stderr,"Error: cannot stat file %s\n",argv[1]);
	  exit(1);
	}
     olen = len = s.st_size;
     Tf = fopen (argv[1],"r");
     if (Tf == NULL)
	{ fprintf (stderr,"Error: cannot open file %s for reading\n",argv[1]);
	  exit(1);
	}
     text = (void*)malloc(len*sizeof(char));
     if (fread(text,1,len,Tf) != len)
	{ fprintf (stderr,"Error: cannot read file %s\n",argv[1]);
	  exit(1);
	}
     fclose(Tf);
     strcpy(fname,argv[1]);
     strcat(fname,".R");
     Rf = fopen (fname,"w");
     if (Rf == NULL)
	{ fprintf (stderr,"Error: cannot open file %s for writing\n",fname);
	  exit(1);
	}
     prepare0 (len);
     len = repair0 (len,Rf);
     if (len == -1)
	{ fprintf (stderr,"Error: cannot write file %s\n",fname);
	  exit(1);
	}
     prepare1 (len);
     if (sC != NULL) len = repair1 (len,Rf);
     if (len == -1)
	{ fprintf (stderr,"Error: cannot write file %s\n",fname);
	  exit(1);
	}
     prepare (len);
     if (repair(Rf) == -1)
	{ fprintf (stderr,"Error: cannot write file %s\n",fname);
	  exit(1);
	}
     if (fclose(Rf) != 0)
	{ fprintf (stderr,"Error: cannot close file %s\n",fname);
	  exit(1);
	}
     strcpy(fname,argv[1]);
     strcat(fname,".C");
     Cf = fopen (fname,"w");
     if (Cf == NULL)
	{ fprintf (stderr,"Error: cannot open file %s for writing\n",fname);
	  exit(1);
	}
     i = 0;
     while (i<u)
        { int cc = (int)C[i];
	  if (fwrite(&cc,sizeof(int),1,Cf) != 1)
             { fprintf (stderr,"Error: cannot write file %s\n",fname);
               exit(1);
             }
          i++; if ((i < u) && (C[i] < 0)) i = -C[i]-1;
        }
     if (fclose(Cf) != 0)
	{ fprintf (stderr,"Error: cannot close file %s\n",fname);
	  exit(1);
	}
     fprintf (stderr,"RePair succeeded\n\n");
     fprintf (stderr,"   Original chars: %lli\n",olen);
     fprintf (stderr,"   Number of rules: %i\n",n-alph);
     fprintf (stderr,"   Final sequence length: %lli\n",c);
     fprintf (stderr,"   Compression ratio: %0.2f%%\n",
			(2.0*(n-alph)+c)*(float)blog(n-1)/(olen*8.0)*100.0);
     exit(0);
   }

