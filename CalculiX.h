/*     CALCULIX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                     */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation; either version 2 of    */
/*     the License, or (at your option) any later version.               */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#define Linux 1
#define IRIX 2
#define IRIX64 3
#define HP 4

#if ARCH == Linux
#define FORTRAN(A,B) A##_  B
#elif ARCH == IRIX || ARCH == IRIX64
#define FORTRAN(A,B) A##_##B
#elif ARCH == HP
#define FORTRAN(A,B) A##B
#endif

#if ARCH == Linux
#define CEE(A,B) A##_  B
#elif ARCH == IRIX || ARCH == IRIX64
#define CEE(A,B) A##_##B
#elif ARCH == HP
#define CEE(A,B) A##B
#endif

#define NNEW(a,b,c) a=(b *)u_calloc((c),sizeof(b),__FILE__,__LINE__,#a)
#define RENEW(a,b,c) a=(b *)u_realloc((b *)(a),(c)*sizeof(b),__FILE__,__LINE__,#a)
#define SFREE(a) u_free(a,__FILE__,__LINE__,#a)

// #define RENEW(a,b,c) a=(b *) realloc((b *)(a),(c)*sizeof(b))

#define DMEMSET(a,b,c,d) for(im=b;im<c;im++)a[im]=d

#ifdef LONGLONG
#define ITG long long
#define ITGFORMAT "lld"
#else
#define ITG int
#define ITGFORMAT "d"
#endif

void FORTRAN(addimdnodecload,(ITG *nodeforc,ITG *i,ITG *imdnode,
             ITG *nmdnode,double *xforc,ITG *ikmpc,ITG *ilmpc,
             ITG *ipompc,ITG *nodempc,ITG *nmpc,ITG *imddof,ITG *nmddof,
             ITG *nactdof,ITG *mi,ITG *imdmpc,ITG *nmdmpc,ITG *imdboun,
	     ITG *nmdboun,ITG *ikboun,ITG *nboun,ITG *ilboun,ITG *ithermal));

void FORTRAN(addimdnodedload,(ITG *nelemload,char *sideload,ITG *ipkon,
             ITG *kon,char *lakon,ITG *i,ITG *imdnode,ITG *nmdnode,
             ITG *ikmpc,ITG *ilmpc,
             ITG *ipompc,ITG *nodempc,ITG *nmpc,ITG *imddof,ITG *nmddof,
             ITG *nactdof,ITG *mi,ITG *imdmpc,ITG *nmdmpc,ITG *imdboun,
	     ITG *nmdboun,ITG *ikboun,ITG *nboun,ITG *ilboun,ITG *ithermal));

void FORTRAN(addizdofcload,(ITG *nodeforc,ITG *ndirforc,ITG *nactdof,
	     ITG *mi,ITG *izdof,ITG *nzdof,ITG *i,ITG *iznode,ITG *nznode,
	     ITG *nk,ITG *imdnode,ITG *nmdnode,double *xforc));

void FORTRAN(addizdofdload,(ITG *nelemload,char *sideload,ITG *ipkon,
             ITG *kon,char *lakon,ITG *nactdof,ITG *izdof,ITG *nzdof,
	     ITG *mi,ITG *i,ITG *iznode,ITG *nznode,ITG *nk,
             ITG *imdnode,ITG *nmdnode));

void FORTRAN(adjustcontactnodes,(char *tieset,ITG *ntie,ITG *itietri,double *cg,
             double *straight,double *co,double *vold,double *xo,double *yo,
             double *zo,double *x,double *y,double *z,ITG *nx,ITG *ny,
             ITG *nz,ITG *istep,ITG *iinc,ITG *iit,ITG *mi,ITG *imastop,
             ITG *nslavnode,ITG *islavnode,char *set,ITG *nset,ITG *istartset,
	     ITG *iendset,ITG *ialset,double *tietol,double *clearini,
	     double *clearslavnode,ITG *itiefac,ITG *ipkon,ITG *kon,
             char *lakon,ITG *islavsurf));

void FORTRAN(allocation,(ITG *nload_,ITG *nforc_,ITG *nboun_,
             ITG *nk_,ITG *ne_,ITG *nmpc_,ITG *nset_,ITG *nalset_,
	     ITG *nmat_,ITG *ntmat_,ITG *npmat_,ITG *norien_,ITG *nam_,
             ITG *nprint_,ITG *mi,ITG *ntrans_,
             char *set,ITG *meminset,
             ITG *rmeminset,ITG *ncs_,ITG *namtot_,ITG *ncmat_,
             ITG *memmpc_,ITG *ne1d,ITG *ne2d,ITG *nflow,
             char *jobnamec,ITG *irstrt,ITG *ithermal,ITG *nener,
             ITG *nstate_,ITG *istep,char *inpc,
             ITG *ipoinp,ITG *inp,ITG *ntie_,ITG *nbody_,
	     ITG *nprop_,ITG *ipoinpc,ITG *nevdamp,ITG *npt_,
	     ITG *nslavsm,ITG *nkon_,ITG *mcs,ITG *mortar,ITG *ifacecount,
	     ITG *nintpoint,ITG *infree,ITG *nheading_,ITG *nobject_));

void FORTRAN(allocont,(ITG *ncont,ITG *ntie,char *tieset,ITG *nset,
             char *set,ITG *istartset,ITG *iendset,ITG *ialset,
	     char *lakon,ITG *ncone,double *tietol,ITG *ismallsliding,
	     char *kind1,char *kind2,ITG *mortar,ITG *istep));

void FORTRAN(applyboun,(ITG *ifaext,ITG *nfaext,ITG *ielfa,ITG *ikboun,
             ITG *ilboun,ITG *nboun,char *typeboun,ITG *nelemload,
             ITG *nload,char *sideload,ITG *isolidsurf,ITG *nsolidsurf,
	     ITG *ifabou,ITG *nfabou,ITG *nface,ITG *nodeboun,ITG *ndirboun,
	     ITG *ikmpc,ITG *ilmpc,char *labmpc,ITG *nmpc,ITG *nactdohinv,
	     ITG *compressible,ITG *iatleastonepressurebc));

void FORTRAN(applympc,(ITG *nface,ITG *ielfa,ITG *is,ITG *ie,ITG *ifabou,
		       ITG *ipompc,double *vfa,double *coefmpc,ITG *nodempc,
                       ITG *ipnei,ITG *neifa,char *labmpc,double *xbounact,
                       ITG *nactdoh,ITG *ifaext,ITG *nfaext));

void FORTRAN(applympc_hfa,(ITG *nface,ITG *ielfa,ITG *is,ITG *ie,ITG *ifabou,
		       ITG *ipompc,double *hfa,double *coefmpc,ITG *nodempc,
                       ITG *ipnei,ITG *neifa,char *labmpc,double *xbounact,
                       ITG *nactdoh));

void arpack(double *co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
	     ITG *ne,
	     ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	     ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
             ITG *nmpc,
	     ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	     ITG *nelemload,char *sideload,double *xload,
	     ITG *nload,
	     ITG *nactdof,
	     ITG *icol,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
	     ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	     ITG *ilboun,
	     double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	     double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,
	     double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
	     ITG **ielorienp,ITG *norien,double *orab,ITG *ntmat_,
	     double *t0,double *t1,double *t1old,
	     ITG *ithermal,double *prestr,ITG *iprestr,
	     double *vold,ITG *iperturb,double *sti,ITG *nzs,  
	     ITG *kode,ITG *mei,double *fei,
	     char *filab,double *eme,
             ITG *iexpl,double *plicon,ITG *nplicon,double *plkcon,
             ITG *nplkcon,
             double **xstatep,ITG *npmat_,char *matname,ITG *mi,
             ITG *ncmat_,ITG *nstate_,double **enerp,char *jobnamec,
             char *output,char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,ITG *nener,ITG *isolver,double *trab,
             ITG *inotr,ITG *ntrans,double *ttime,double *fmpc,
	     char *cbody,ITG *ibody,double *xbody,ITG *nbody,double *thicke,
	     ITG *nslavs,double *tietol,ITG *nkon,ITG *mpcinfo,ITG *ntie,
	     ITG *istep,ITG *mcs,ITG *ics,char *tieset,
             double *cs,ITG *nintpoint,ITG *mortar,ITG *ifacecount,
	     ITG **islavsurfp,double **pslavsurfp,double **clearinip,
	     ITG *nmat,char *typeboun,ITG *ielprop,double *prop);

void arpackbu(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	     ITG *ne,
	     ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	     ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
             ITG *nmpc,
	     ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	     ITG *nelemload,char *sideload,double *xload,
	     ITG *nload,
	     ITG *nactdof,
	     ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
	     ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	     ITG *ilboun,
	     double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	     double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	     ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	     double *t0,double *t1,double *t1old,
	     ITG *ithermal,double *prestr,ITG *iprestr,
	     double *vold,ITG *iperturb,double *sti,ITG *nzs,
	     ITG *kode,ITG *mei,double *fei,
             char *filab,double *eme,
             ITG *iexpl,double *plicon,ITG *nplicon,double *plkcon,
             ITG *nplkcon,
             double *xstate,ITG *npmat_,char *matname,ITG *mi,
             ITG *ncmat_,ITG *nstate_,double *ener,char *output,
             char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,ITG *nener,ITG *isolver,double *trab,
             ITG *inotr,ITG *ntrans,double *ttime,double *fmpc,
	     char *cbody,ITG *ibody,double *xbody,ITG *nbody,
	     double *thicke,char *jobnamec,ITG *nmat,ITG *ielprop,
             double *prop);

void arpackcs(double *co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
	     ITG *ne,
	     ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	     ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
             ITG *nmpc,
	     ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	     ITG *nelemload,char *sideload,double *xload,
	     ITG *nload,ITG *nactdof,
	     ITG *icol,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
	     ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	     ITG *ilboun,
	     double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	     double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
	     ITG **ielorienp,ITG *norien,double *orab,ITG *ntmat_,
	     double *t0,double *t1,double *t1old,
	     ITG *ithermal,double *prestr,ITG *iprestr,
	     double *vold,ITG *iperturb,double *sti,ITG *nzs, 
	     ITG *kode,ITG *mei,double *fei,
	     char *filab,double *eme,
             ITG *iexpl,double *plicon,ITG *nplicon,double *plkcon,
             ITG *nplkcon,
             double **xstatep,ITG *npmat_,char *matname,ITG *mi,
             ITG *ics,double *cs,ITG *mpcend,ITG *ncmat_,ITG *nstate_,
             ITG *mcs,ITG *nkon,double **enerp,char *jobnamec,
             char *output,char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,ITG *nener,ITG *isolver,double *trab,
             ITG *inotr,ITG *ntrans,double *ttime,double *fmpc,
	     char *cbody,ITG *ibody,double *xbody,ITG *nbody,
             ITG *nevtot,double *thicke,ITG *nslavs,double *tietol,
	     ITG *mpcinfo,ITG *ntie,ITG *istep,
	     char *tieset,ITG *nintpoint,ITG *mortar,ITG *ifacecount,
	     ITG **islavsurfp,double **pslavsurfp,double **clearinip,
	     ITG *nmat,char *typeboun,ITG *ielprop,double *prop);

void FORTRAN(assigndomtonodes,(ITG *ne,char *lakon,ITG *ipkon,ITG *kon,
             ITG *ielmat,ITG *inomat,double *elcon,ITG *ncmat_,ITG *ntmat_,
	     ITG *mi,ITG *ne2));

void FORTRAN(basis,(double *x,double *y,double *z,double *xo,double *yo,
                    double *zo,ITG *nx,ITG *ny,ITG *nz,double *planfa,
                    ITG *ifatet,ITG *nktet,ITG *netet,double *field,
                    ITG *nfield,double *cotet,ITG *kontyp,ITG *ipkon,
                    ITG *kon,ITG *iparent,double *xp,double *yp,double *zp,
                    double *value,double *ratio,ITG *iselect,ITG *nselect,
                    ITG *istartset,ITG *iendset,ITG *ialset,ITG *imastset,
                    ITG *ielemnr,ITG *nterms,ITG *konl));

void biosav(ITG *ipkon,ITG *kon,char *lakon,ITG *ne,double *co,
	    double *qfx,double *h0,ITG *mi,ITG *inomat,ITG *nk);

void FORTRAN(biotsavart,(ITG *ipkon,ITG *kon,char *lakon,ITG *ne,double *co,
                         double *qfx,double *h0,ITG *mi,ITG *nka,ITG *nkb));

void *biotsavartmt(ITG *i);

void FORTRAN(bodyforce,(char *cbody,ITG *ibody,ITG *ipobody,ITG *nbody,
             char *set,ITG *istartset,ITG *iendset,ITG *ialset,
             ITG *inewton,ITG *nset,ITG *ifreebody,ITG *k));

void FORTRAN(calcbody,(ITG *nef,double *body,ITG *ipobody,ITG *ibody,
		       double *xbody,double *coel,double *vel,char *lakon,
                       ITG *nactdohinv));

void FORTRAN(calcguesstincf,(ITG *nface,double *dmin,double *vfa,
                             double *umfa,double *cvfa,double *hcfa,
			     ITG *ithermal,double *tincfguess,
                             ITG *compressible));

void FORTRAN(calcinitialflux,(double *area,double *vfa,double *xxn,
                              ITG *ipnei,ITG *nef,ITG *neifa,char *lakonf,
                              double *flux));

void FORTRAN(calccvel,(ITG *nef,double *vel,double *shcon,
		       ITG *nshcon,ITG *ielmatf,ITG *ntmat_,
                       ITG *mi,double *cvel,double *physcon));

void FORTRAN(calccvelcomp,(ITG *nef,double *vel,double *shcon,
		       ITG *nshcon,ITG *ielmatf,ITG *ntmat_,
                       ITG *mi,double *cvel,double *physcon));

void FORTRAN(calccvfa,(ITG *nface,double *vfa,double *shcon,
			ITG *nshcon,ITG *ielmatf,ITG *ntmat_,
		       ITG *mi,ITG *ielfa,double *cvfa,double *physcon));

void FORTRAN(calccvfacomp,(ITG *nface,double *vfa,double *shcon,
			ITG *nshcon,ITG *ielmatf,ITG *ntmat_,
		       ITG *mi,ITG *ielfa,double *cvfa,double *physcon));

void FORTRAN(calcgamma,(ITG *nface,ITG *ielfa,double *vel,double *gradvel,
                        double *gamma,double *xlet,double *xxn,double *xxj,
                        ITG *ipnei,double *betam,ITG *nef,double *flux));

void FORTRAN(calcgammat,(ITG *nface,ITG *ielfa,double *vel,double *gradtfa,
                        double *gammat,double *xlet,double *xxn,double *xxj,
			ITG *ipnei,double *betam,ITG *nef));

void FORTRAN(calcgradpel,(ITG *ne,char *lakon,ITG *ipnei,double *vfa,
                          double *area,double *xxn,double *gradpel,ITG *neifa,
			  double *volume));

void FORTRAN(calcgradtel,(ITG *ne,char *lakon,ITG *ipnei,double *vfa,
                          double *area,double *xxn,double *gradtel,ITG *neifa,
			  double *volume));

void FORTRAN(calcgradvel,(ITG *ne,char *lakon,ITG *ipnei,double *vfa,
			double *area,double *xxn,double *gradv,ITG *neifa,
                        double *volume));

void FORTRAN(calchcel,(double *vel,double *cocon,ITG *ncocon,ITG *ielmatf,
                       ITG *ntmat_,ITG *mi,double *hcel,ITG *nef));

void FORTRAN(calchcfa,(ITG *nface,double *vfa,double *cocon,
			ITG *ncocon,ITG *ielmatf,ITG *ntmat_,
                        ITG *mi,ITG *ielfa,double *hcfa));

void FORTRAN(calch0interface,(ITG *nmpc,ITG *ipompc,ITG *nodempc,
                              double *coefmpc,double *h0));
		      
void FORTRAN(calcmac,(ITG *neq,double *z,double *zz,ITG *nev,double *mac,
		      double* maccpx,ITG *istartnmd,ITG *iendnmd,ITG *nmd,
		      ITG *cyclicsymmetry,ITG *neqact,double *bett,
		      double *betm));

void FORTRAN(calcmass,(ITG *ipkon,char *lakon,ITG *kon,double *co,ITG *mi,
             ITG *nelem,ITG *ne,double *thicke,ITG *ielmat,
             ITG *nope,double *t0,double *t1,double *rhcon,
	     ITG *nrhcon,ITG *ntmat_,ITG *ithermal,double *csmass,
	     ITG *ielprop,double *prop));

void FORTRAN(calcmatwavspeed,(ITG *ne0,double *elcon,ITG *nelcon,
            double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,double *orab,
            ITG *ntmat_,ITG *ithermal,double *alzero,double *plicon,
            ITG *nplicon,double *plkcon,ITG *nplkcon,ITG *npmat_,ITG *mi,
            double *dtime,double *xstiff,ITG *ncmat_,double *vold,ITG *ielmat,
            double *t0,double *t1,char *matname,char *lakon,
	    double *xmatwavespeed,ITG *nmat,ITG *ipkon));

void FORTRAN(calcpel,(ITG *ne,ITG *nactdoh,double *vel,double *b,ITG *nef));

void calcresidual(ITG *nmethod,ITG *neq,double *b,double *fext,double *f,
        ITG *iexpl,ITG *nactdof,double *aux2,double *vold,
        double *vini,double *dtime,double *accold,ITG *nk,double *adb,
        double *aub,ITG *icol,ITG *irow,ITG *nzl,double *alpha,
	double *fextini,double *fini,ITG *islavnode,ITG *nslavnode,
        ITG *mortar,ITG *ntie,
        double *f_cm,double *f_cs,ITG *mi,ITG *nzs,ITG *nasym,
	ITG *idamping,double *veold,double *adc,double *auc,double *cvini,
        double *cv);

void calcresidual_em(ITG *nmethod,ITG *neq,double *b,double *fext,double *f,
        ITG *iexpl,ITG *nactdof,double *aux1,double *aux2,double *vold,
        double *vini,double *dtime,double *accold,ITG *nk,double *adb,
        double *aub,ITG *icol,ITG *irow,ITG *nzl,double *alpha,
	double *fextini,double *fini,ITG *islavnode,ITG *nslavnode,
        ITG *mortar,ITG *ntie,
	double *f_cm,double *f_cs,ITG *mi,ITG *nzs,ITG *nasym,ITG *ithermal);

void FORTRAN(calcrhoel,(ITG *nef,double *vel,double *rhcon,
			ITG *nrhcon,ITG *ielmatf,ITG *ntmat_,ITG *ithermal,
                        ITG *mi));

void FORTRAN(calcrhoelcomp,(ITG *nef,double *vel,double *shcon,
			ITG *ielmatf,ITG *ntmat_,ITG *mi));

void FORTRAN(calcrhofa,(ITG *nface,double *vfa,double *rhcon,
			ITG *nrhcon,ITG *ielmatf,ITG *ntmat_,ITG *ithermal,
                        ITG *mi,ITG *ielfa));

void FORTRAN(calcrhofacomp,(ITG *nface,double *vfa,double *shcon,
			ITG *ielmatf,ITG *ntmat_,ITG *mi,ITG *ielfa));

void FORTRAN(calcstabletimeinccont,(ITG *ne,char *lakon,ITG *kon,ITG *ipkon,
             ITG *mi,ITG *ielmat,double *elcon,ITG *mortar,double *adb,
             double *alpha,ITG *nactdof,double *springarea,ITG *ne0,
             ITG *ntmat_,ITG *ncmat_,double *dtcont));

void FORTRAN(calcstabletimeincvol,(ITG *ne0,char *lakon,double *co,
	     ITG *kon,ITG *ipkon,ITG *mi,ITG *ielmat,double *dtvol,
             double *alpha,double *xmatwavespeed));

void FORTRAN(calcstressheatflux,(double *sti,double *umel,double *gradvel,
             double *qfx,double *hcel,double *gradtel,ITG *nef,ITG *isti,
	     ITG *iqfx,ITG *mi));

void FORTRAN(calctel,(ITG *ne,ITG *nactdoh,double *vel,double *b,ITG *nef));

void FORTRAN(calcumel,(ITG *nef,double *vel,double *shcon,
		       ITG *nshcon,ITG *ielmatf,ITG *ntmat_,ITG *ithermal,
                       ITG *mi,double *umel));

void FORTRAN(calcumfa,(ITG *nface,double *vfa,double *shcon,
			ITG *nshcon,ITG *ielmatf,ITG *ntmat_,ITG *ithermal,
                        ITG *mi,ITG *ielfa,double *umfa));

void FORTRAN(calcvel,(ITG *ne,ITG *nactdoh,double *vel,double *b,
                      ITG *neq,ITG *nef));

void FORTRAN(calcview,(char *sideload,double *vold,double *co,
             double *pmid,double *e1,double *e2,double *e3,
	     ITG *kontri,ITG *nloadtr,double *adview,double *auview,
             double *dist,ITG *idist,double *area,ITG *ntrit,ITG *mi,ITG *jqrad,
	     ITG *irowrad,ITG *nzsrad,double *sidemean,ITG *ntria,
             ITG *ntrib));

void *calcviewmt(ITG *i);

void FORTRAN(calinput,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *nkon,ITG *ne,
	       ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	       ITG *ipompc,ITG *nodempc,double *coefmpc,ITG *nmpc,
	       ITG *nmpc_,ITG *nodeforc,ITG *ndirforc,double *xforc,
	       ITG *nforc,ITG *nforc_,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,ITG *nload_,
	       ITG *nprint,char *prlab,char *prset,ITG *mpcfree,ITG *nboun_,
	       ITG *mei,char *set,ITG *istartset,
	       ITG *iendset,ITG *ialset,ITG *nset,ITG *nalset,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,double *t0,
	       double *t1,char *matname,
	       ITG *ielmat,char *orname,double *orab,ITG *ielorien,
	       char *amname,double *amta,ITG *namta,ITG *nam,
	       ITG *nmethod,ITG *iamforc,ITG *iamload,
	       ITG *iamt1,ITG *ithermal,ITG *iperturb,
	       ITG *istat,ITG *istep,ITG *nmat,
	       ITG *ntmat_,ITG *norien,double *prestr,ITG *iprestr,
	       ITG *isolver,double *fei,double *veold,double *timepar,
	       double *xmodal,char *filab,ITG *jout,ITG *nlabel,
	       ITG *idrct,ITG *jmax,ITG *iexpl,double *alpha,ITG *iamboun,
	       double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
	       ITG *iplas,ITG *npmat_,ITG *mi,ITG *nk_,
	       double *trab,ITG *inotr,ITG *ntrans,ITG *ikboun,
               ITG *ilboun,ITG *ikmpc,ITG *ilmpc,ITG *ics,
	       double *dcs,ITG *ncs_,ITG *namtot_,double *cs,
               ITG *nstate_,ITG *ncmat_,ITG *iumat,ITG *mcs,
               char *labmpc,ITG *iponor,double *xnor,ITG *knor,
	       double *thickn,double *thicke,ITG *ikforc,ITG *ilforc,
               double *offset,ITG *iponoel,ITG *inoel,ITG *rig,
               ITG *infree,ITG *nshcon,double *shcon,double *cocon,
               ITG *ncocon,double *physcon,ITG *nflow,double *ctrl,
               ITG *maxlenmpc,ITG *ne1d,ITG *ne2d,ITG *nener,
               double *vold,ITG *nodebounold,
               ITG *ndirbounold,double *xbounold,double *xforcold,
               double *xloadold,double *t1old,double *eme,
               double *sti,double *ener,
               double *xstate,char *jobnamec,ITG *irstrt,
               double *ttime,double *qaold,
               char *output,char *typeboun,char *inpc,
               ITG *ipoinp,ITG *inp,char *tieset,double *tietol,
               ITG *ntie,double *fmpc,char *cbody,ITG *ibody,double *xbody,
               ITG *nbody,ITG *nbody_,double *xbodyold,ITG *nam_,
               ITG *ielprop,ITG *nprop,ITG *nprop_,double *prop,
	       ITG *itpamp,ITG *iviewfile,ITG *ipoinpc,ITG *cfd,
	       ITG *nslavs,double *t0g,double *t1g,ITG *network,
	       ITG *cyclicsymmetry,ITG *idefforc,ITG *idefload,
               ITG *idefbody,ITG *mortar,ITG *ifacecount,ITG *islavsurf,
	       double *pslavsurf,double *clearini,char *heading,
               ITG *iaxial,ITG *nobject,char *objectset));    

void cascade(ITG *ipompc,double **coefmpcp,ITG **nodempcp,ITG *nmpc,
   ITG *mpcfree,ITG *nodeboun,ITG *ndirboun,ITG*nboun,ITG*ikmpc,
   ITG *ilmpc,ITG *ikboun,ITG *ilboun,ITG *mpcend,ITG *mpcmult,
   char *labmpc,ITG *nk,ITG *memmpc_,ITG *icascade,ITG *maxlenmpc,
   ITG *callfrommain,ITG *iperturb,ITG *ithermal);

void FORTRAN(cataloguenodes,(ITG *iponofa,ITG *inofa,ITG *ifreefa,ITG *ielfa,
             ITG *ifaboun,ITG *ipkon,ITG *kon,char *lakon,ITG *nface,ITG *ne));

ITG cgsolver(double *A,double *x,double *b,ITG neq,ITG len,ITG *ia,ITG *iz,
				double *eps,ITG *niter,ITG precFlg);

void checkconvergence(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	  ITG *ne,double *stn,ITG *nmethod,
	  ITG *kode,char *filab,double *een,double *t1act,
          double *time,double *epn,ITG *ielmat,char *matname,
          double *enern,double *xstaten,ITG *nstate_,ITG *istep,
          ITG *iinc,ITG *iperturb,double *ener,ITG *mi,char *output,
          ITG *ithermal,double *qfn,ITG *mode,ITG *noddiam,double *trab,
          ITG *inotr,ITG *ntrans,double *orab,ITG *ielorien,ITG *norien,
          char *description,double *sti,
	  ITG *icutb,ITG *iit,double *dtime,double *qa,double *vold,
          double *qam,double *ram1,double *ram2,double *ram,
          double *cam,double *uam,ITG *ntg,double *ttime,
          ITG *icntrl,double *theta,double *dtheta,double *veold,
          double *vini,ITG *idrct,double *tper,ITG *istab,double *tmax,
	  ITG *nactdof,double *b,double *tmin,double *ctrl,double *amta,
          ITG *namta,ITG *itpamp,ITG *inext,double *dthetaref,ITG *itp,
          ITG *jprint,ITG *jout,ITG *uncoupled,double *t1,ITG *iitterm,
          ITG *nelemload,ITG *nload,ITG *nodeboun,ITG *nboun,ITG *itg,
	  ITG *ndirboun,double *deltmx,ITG *iflagact,char *set,ITG *nset,
	  ITG *istartset,ITG *iendset,ITG *ialset,double *emn,double *thicke,
	  char *jobnamec,ITG *mortar,ITG *nmat,ITG *ielprop,double *prop,
	  ITG *ialeatoric,ITG *kscale,
          double *energy,double *allwk,double *energyref,
          double *emax,double *enres,double *enetoll,double *energyini,
          double *allwkini,double *temax,double *reswk,ITG *ne0,
          ITG *neini,double *dampwk,double *dampwkini,double *energystartstep);

void checkconvnet(ITG *icutb,ITG *iin,
		  double *qamt,double *qamf,double *qamp,
		  double *ram1t,double *ram1f,double *ram1p,
		  double *ram2t,double *ram2f,double *ram2p,
		  double *ramt,double *ramf,double *ramp,
		  ITG *icntrl,double *dtheta,double *ctrl,
                  double *uama,double *ram1a,double *ram2a,double *rama,
                  double *vamt,double *vamf,double *vamp,double *vama,
                  double *qa);

void checkinclength(double *time,double *ttime,double *theta,double *dtheta,
          ITG *idrct,double *tper,double *tmax,double *tmin,double *ctrl,
          double *amta,ITG *namta,ITG *itpamp,ITG *inext,double *dthetaref,
	  ITG *itp,ITG *jprint,ITG *jout);
         
void FORTRAN(checkimpacts,(ITG *ne,ITG *neini,double *temax,
                   double *sizemaxinc,double *energyref,double *tmin,
                   double *tper,
                   ITG *idivergence,ITG *idirinctime,ITG *istab,
                   double *dtheta,double *enres,double *energy,
		   double *energyini,double *allwk,double *allwkini,
		   double *dampwk,double *dampwkini,double *emax,
		   ITG *mortar,double *maxdecay,double *enetoll));

void FORTRAN(checktime,(ITG *itpamp,ITG *namta,double *tinc,double *ttime,
	     double *amta,double *tmin,ITG *inext,ITG *itp,ITG *istep));

void FORTRAN(checktruecontact,(ITG *ntie,char *tieset,double *tietol,
	     double *elcon,ITG *itruecontact,ITG *ncmat_,ITG *ntmat_));

void FORTRAN(closefile,());

void FORTRAN(closefilefluid,());

void compfluid(double **cop,ITG *nk,ITG **ipkonp,ITG **konp,char **lakonp,
    char **sideface,ITG *ifreestream,
    ITG *nfreestream,ITG *isolidsurf,ITG *neighsolidsurf,
    ITG *nsolidsurf,ITG **iponoel,ITG **inoel,ITG *nshcon,double *shcon,
    ITG *nrhcon,double *rhcon,double **voldp,ITG *ntmat_,ITG *nodeboun,
    ITG *ndirboun,ITG *nboun,ITG **ipompcp,ITG **nodempcp,ITG *nmpc,
    ITG **ikmpcp,ITG **ilmpcp,ITG *ithermal,ITG *ikboun,ITG *ilboun,
    ITG *turbulent,ITG *isolver,ITG *iexpl,double *vcontu,double *ttime,
    double *time,double *dtime,ITG *nodeforc,ITG *ndirforc,double *xforc,
    ITG *nforc,ITG *nelemload,char *sideload,double *xload,ITG *nload,
    double *xbody,ITG *ipobody,ITG *nbody,ITG *ielmatf,char *matname,
    ITG *mi,ITG *ncmat_,double *physcon,ITG *istep,ITG *iinc,
    ITG *ibody,double *xloadold,double *xboun,
    double **coefmpcp,ITG *nmethod,double *xforcold,double *xforcact,
    ITG *iamforc,ITG *iamload,double *xbodyold,double *xbodyact,
    double *t1old,double *t1,double *t1act,ITG *iamt1,double *amta,
    ITG *namta,ITG *nam,double *ampli,double *xbounold,double *xbounact,
    ITG *iamboun,ITG *itg,ITG *ntg,char *amname,double *t0,ITG **nelemface,
    ITG *nface,double *cocon,ITG *ncocon,double *xloadact,double *tper,
    ITG *jmax,ITG *jout,char *set,ITG *nset,ITG *istartset,
    ITG *iendset,ITG *ialset,char *prset,char *prlab,ITG *nprint,
    double *trab,ITG *inotr,ITG *ntrans,char *filab,char **labmpcp,
    double *sti,ITG *norien,double *orab,char *jobnamef,char *tieset,
    ITG *ntie,ITG *mcs,ITG *ics,double *cs,ITG *nkon,ITG *mpcfree,
    ITG *memmpc_,double **fmpcp,ITG *nef,ITG **inomat,double *qfx,
    ITG *neifa,ITG *neiel,ITG *ielfa,ITG *ifaext,double *vfa,double *vel,
    ITG *ipnei,ITG *nflnei,ITG *nfaext,char *typeboun,ITG *neij,
    double *tincf,ITG *nactdoh,ITG *nactdohinv,ITG *ielorien,
    char *jobnamec);

void FORTRAN(complete_hel,(ITG *neq,double *b,double *hel,double *ad,
             double *au,ITG *jq,ITG *irow,ITG *nzs));

void complexfreq(double **cop,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,ITG *ne,
	       ITG **nodebounp,ITG **ndirbounp,double **xbounp,ITG *nboun,
	       ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,double *xforc,
               ITG *nforc,ITG *nelemload,char *sideload,double *xload,
	       ITG *nload,
	       ITG **nactdofp,ITG *neq,ITG *nzl,ITG *icol,ITG *irow,
	       ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,ITG **ikbounp,
	       ITG **ilbounp,double *elcon,ITG *nelcon,double *rhcon,
	       ITG *nrhcon,double *cocon,ITG *ncocon,
               double *alcon,ITG *nalcon,double *alzero,
               ITG **ielmatp,ITG **ielorienp,ITG *norien,double *orab,
               ITG *ntmat_,double **t0p,
	       double **t1p,ITG *ithermal,double *prestr,ITG *iprestr,
	       double **voldp,ITG *iperturb,double **stip,ITG *nzs,
	       double *timepar,double *xmodal,
	       double **veoldp,char *amname,double *amta,
	       ITG *namta,ITG *nam,ITG *iamforc,ITG *iamload,
	       ITG **iamt1p,ITG *jout,
	       ITG *kode,char *filab,double **emep,double *xforcold,
	       double *xloadold,
               double **t1oldp,ITG **iambounp,double **xbounoldp,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstate,ITG *npmat_,char *matname,ITG *mi,
               ITG *ncmat_,ITG *nstate_,double **enerp,char *jobnamec,
               double *ttime,char *set,ITG *nset,ITG *istartset,
               ITG *iendset,ITG **ialsetp,ITG *nprint,char *prlab,
               char *prset,ITG *nener,double *trab,
               ITG **inotrp,ITG *ntrans,double **fmpcp,char *cbody,ITG *ibody,
               double *xbody,ITG *nbody,double *xbodyold,ITG *istep,
               ITG *isolver,ITG *jq,char *output,ITG *mcs,ITG *nkon,
               ITG *mpcend,ITG *ics,double *cs,ITG *ntie,char *tieset,
               ITG *idrct,ITG *jmax,
	       double *ctrl,ITG *itpamp,double *tietol,ITG *nalset,
	       ITG *ikforc,ITG *ilforc,double *thicke,
	       char *jobnamef,ITG *mei,ITG *nmat,ITG *ielprop,double *prop);

void contact(ITG *ncont,ITG *ntie,char *tieset,ITG *nset,char *set,
	     ITG *istartset,ITG *iendset,ITG *ialset,ITG *itietri,
	     char *lakon,ITG *ipkon,ITG *kon,ITG *koncont,ITG *ne,
	     double *cg,double *straight,ITG *ifree,double *co,
	     double *vold,ITG *ielmat,double *cs,double *elcon,
             ITG *istep,ITG *iinc,ITG *iit,ITG *ncmat_,ITG *ntmat_,
             ITG *ne0,double *vini,
             ITG *nmethod,ITG *nmpc,ITG *mpcfree,ITG *memmpc_,
             ITG **ipompcp,char **labmpcp,ITG **ikmpcp,ITG **ilmpcp,
             double **fmpcp,ITG **nodempcp,double **coefmpcp,
             ITG *iperturb,ITG *ikboun,ITG *nboun,ITG *mi,ITG *imastop,
             ITG *nslavnode,ITG *islavnode,ITG *islavsurf,ITG *itiefac,
             double *areaslav,ITG *iponoels,ITG *inoels,double *springarea,
             double *tietol,double *reltime,ITG *imastnode,ITG *nmastnode,
             double *xmastnor,char *filab,ITG *mcs,
             ITG *ics,ITG *nasym,double *xnoels,ITG *mortar,
             double *pslavsurf,double *pmastsurf,double *clearini,
             double *theta,double *xstateini,double *xstate,ITG *nstate_,
             ITG *icutb,ITG *ialeatoric);

void FORTRAN(convert2slapcol,(double *au,double *ad,ITG *irow,ITG *ia,
                              ITG *jq,ITG *ja,ITG *nzs,ITG *neq,double *aua));
    
void FORTRAN(coriolissolve,(double *cc,ITG *nev,double *aa,double *bb,
             double *xx,double *eiga,double *eigb,double *eigxx,
	     ITG *iter,double *d,double *temp));

void FORTRAN(correctvel,(double *hel,double *adv,double *vfa,ITG *ipnei,
                         double *area,double *bv,double *xxn,ITG *neifa,
			 char *lakon,ITG *ne,ITG *neq));

void FORTRAN(correctvfa,(ITG *nface,ITG *ielfa,double *area,double *vfa,
		       double *ap,double *bp,
                       double *xxn,ITG *ifabou,ITG *ipnei,
                       ITG *nef,ITG *neifa,double *hfa,double *vel,
		       double *xboun,char *lakonf,double *flux));

void FORTRAN(createinterfacempcs,(ITG *imastnode,double *xmastnor,
	     ITG *nmastnode,ITG *ikmpc,ITG *ilmpc,ITG *nmpc,ITG *ipompc,
             ITG *nodempc,double *coefmpc,char *labmpc,ITG *mpcfree,
             ITG *ikboun,ITG *nboun));

void FORTRAN(createinum,(ITG *ipkon,ITG *inum,ITG *kon,char *lakon,ITG *nk,
             ITG *ne,char *cflag,ITG *nelemload,ITG *nload,ITG *nodeboun,
             ITG *nboun,ITG *ndirboun,ITG *ithermal,double *co,
	     double *vold,ITG *mi,ITG *ielmat));

void FORTRAN(createmddof,(ITG *imddof,ITG *nmddof,ITG *istartset,
       ITG *iendset,ITG *ialset,ITG *nactdof,ITG *ithermal,ITG *mi,
       ITG *imdnode,ITG *nmdnode,ITG *ikmpc,ITG *ilmpc,ITG *ipompc,
       ITG *nodempc,ITG *nmpc,ITG *imdmpc,
       ITG *nmdmpc,ITG *imdboun,ITG *nmdboun,ITG *ikboun,ITG *nboun,
       ITG *nset,ITG *ntie,char *tieset,char *set,char *lakon,ITG *kon,
       ITG *ipkon,char *labmpc,ITG *ilboun,char *filab,char *prlab,
       char *prset,ITG *nprint,ITG *ne,ITG *cyclicsymmetry));

void FORTRAN(createmdelem,(ITG *imdnode,ITG *nmdnode,double *xforc,
             ITG *ikmpc,ITG *ilmpc,ITG *ipompc,ITG *nodempc,ITG *nmpc,
             ITG *imddof,ITG *nmddof,ITG *nactdof,ITG *mi,ITG *imdmpc,
             ITG *nmdmpc,ITG *imdboun,ITG *nmdboun,ITG *ikboun,ITG *nboun,
             ITG *ilboun,ITG *ithermal,ITG *imdelem,ITG *nmdelem,
             ITG *iponoel,ITG *inoel,char *prlab,char *prset,ITG *nprint,
             char *lakon,char *set,ITG *nset,ITG *ialset,ITG *ipkon,
             ITG *kon,ITG *istartset,ITG *iendset,ITG *nforc,
             ITG *ikforc,ITG *ilforc));

void FORTRAN(createtiedsurfs,(ITG *nodface,ITG *ipoface,char *set,
             ITG *istartset,ITG *iendset,ITG *ialset,char *tieset,
             ITG *inomat,ITG *ne,ITG *ipkon,char *lakon,ITG *kon,
	     ITG *ntie,double *tietol,ITG *nalset,ITG *nk,ITG *nset,
             ITG *iactive));

void FORTRAN(dattime,(char *date,char *clock));

void CEE(ddotc,(ITG *n,double *dx,ITG *incx,double *dy,ITG *incy,
                double *funcddot));

void *ddotc1mt(ITG *i);

void dfdbj(double *bcont,double **dbcontp,ITG *neq,ITG *nope,
	   ITG *konl,ITG *nactdof,double *s,double *z,ITG *ikmpc,
	   ITG *ilmpc,ITG *ipompc,ITG *nodempc,ITG *nmpc,
	   double *coefmpc,double *fnl,ITG *nev,
	   ITG **ikactcontp,ITG **ilactcontp,ITG *nactcont,ITG *nactcont_,
           ITG *mi,ITG *cyclicsymmetry,ITG *izdof,ITG *nzdof);
      
void FORTRAN(dgesv,(ITG *nteq,ITG *nhrs,double *ac,ITG *lda,ITG *ipiv,
                     double *bc,ITG *ldb,ITG *info)); 

void FORTRAN(dgetrs,(char *trans,ITG *nteq,ITG *nrhs,double *ac,ITG *lda,
		      ITG *ipiv,double *bc,ITG *ldb,ITG *info));

void FORTRAN(drfftf,(ITG *ndata,double *r,double *wsave,ITG *isave));

void FORTRAN(drffti,(ITG *ndata,double *wsave,ITG *isave));

void FORTRAN(dnaupd,(ITG *ido,char *bmat,ITG *n,char *which,ITG *nev,
	     double *tol,double *resid,ITG *ncv,double *z,ITG *ldz,
	     ITG *iparam,ITG *ipntr,double *workd,double *workl,
	     ITG *lworkl,ITG *info));

void FORTRAN(dsaupd,(ITG *ido,char *bmat,ITG *n,char *which,ITG *nev,
	     double *tol,double *resid,ITG *ncv,double *z,ITG *ldz,
	     ITG *iparam,ITG *ipntr,double *workd,double *workl,
	     ITG *lworkl,ITG *info));

void FORTRAN(dneupd,(ITG *rvec,char *howmny,ITG *select,double *d,
	     double *di,double *z,ITG *ldz,double *sigma,double *sigmai,
             double *workev,char *bmat,ITG *neq,char *which,
	     ITG *nev,double *tol,double *resid,ITG *ncv,double *v,
	     ITG *ldv,ITG *iparam,ITG *ipntr,double *workd,
	     double *workl,ITG *lworkl,ITG *info));

void FORTRAN(dseupd,(ITG *rvec,char *howmny,ITG *select,double *d,double *z,
	     ITG *ldz,double *sigma,char *bmat,ITG *neq,char *which,
	     ITG *nev,double *tol,double *resid,ITG *ncv,double *v,
	     ITG *ldv,ITG *iparam,ITG *ipntr,double *workd,
	     double *workl,ITG *lworkl,ITG *info));

void FORTRAN(dsort,(double *dx,ITG *iy,ITG *n,ITG *kflag));

void dyna(double **cop,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,ITG *ne,
	       ITG **nodebounp,ITG **ndirbounp,double **xbounp,ITG *nboun,
	       ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,double *xforc,
               ITG *nforc,ITG *nelemload,char *sideload,double *xload,
	       ITG *nload,
	       ITG **nactdofp,ITG *neq,ITG *nzl,ITG *icol,ITG *irow,
	       ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,ITG **ikbounp,
	       ITG **ilbounp,double *elcon,ITG *nelcon,double *rhcon,
	       ITG *nrhcon,double *cocon,ITG *ncocon,
               double *alcon,ITG *nalcon,double *alzero,
               ITG **ielmatp,ITG **ielorienp,ITG *norien,double *orab,
               ITG *ntmat_,double **t0p,
	       double **t1p,ITG *ithermal,double *prestr,ITG *iprestr,
	       double **voldp,ITG *iperturb,double **stip,ITG *nzs,
	       double *timepar,double *xmodal,
	       double **veoldp,char *amname,double *amta,
	       ITG *namta,ITG *nam,ITG *iamforc,ITG *iamload,
	       ITG **iamt1p,ITG *jout,
	       ITG *kode,char *filab,double **emep,double *xforcold,
	       double *xloadold,
               double **t1oldp,ITG **iambounp,double **xbounoldp,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double **xstatep,ITG *npmat_,char *matname,ITG *mi,
               ITG *ncmat_,ITG *nstate_,double **enerp,char *jobnamec,
               double *ttime,char *set,ITG *nset,ITG *istartset,
               ITG *iendset,ITG **ialsetp,ITG *nprint,char *prlab,
               char *prset,ITG *nener,double *trab,
               ITG **inotrp,ITG *ntrans,double **fmpcp,char *cbody,ITG *ibody,
               double *xbody,ITG *nbody,double *xbodyold,ITG *istep,
               ITG *isolver,ITG *jq,char *output,ITG *mcs,ITG *nkon,
               ITG *mpcend,ITG *ics,double *cs,ITG *ntie,char *tieset,
               ITG *idrct,ITG *jmax,
	       double *ctrl,ITG *itpamp,double *tietol,ITG *nalset,
	       ITG *ikforc,ITG *ilforc,double *thicke,
	       ITG *nslavs,ITG *nmat,char *typeboun,ITG *ielprop,double *prop);

void dynacont(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
	      ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	      ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
	      ITG *nmpc,ITG *nodeforc,ITG *ndirforc,double *xforc,
	      ITG *nforc,ITG *nelemload,char *sideload,double *xload,
	      ITG *nload,
	      ITG *nactdof,ITG *neq,ITG *nzl,ITG *icol,ITG *irow,
	      ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	      ITG *ilboun,double *elcon,ITG *nelcon,double *rhcon,
	      ITG *nrhcon,double *cocon,ITG *ncocon,
	      double *alcon,ITG *nalcon,double *alzero,
	      ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,
	      ITG *ntmat_,double *t0,
	      double *t1,ITG *ithermal,double *prestr,ITG *iprestr,
	      double *vold,ITG *iperturb,double *sti,ITG *nzs,
	      double *tinc,double *tper,double *xmodal,
	      double *veold,char *amname,double *amta,
	      ITG *namta,ITG *nam,ITG *iamforc,ITG *iamload,
	      ITG *iamt1,ITG *jout,char *filab,double *eme,double *xforcold,
	      double *xloadold,
	      double *t1old,ITG *iamboun,double *xbounold,ITG *iexpl,
	      double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
	      double *xstate,ITG *npmat_,char *matname,ITG *mi,
	      ITG *ncmat_,ITG *nstate_,double *ener,char *jobnamec,
	      double *ttime,char *set,ITG *nset,ITG *istartset,
	      ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
	      char *prset,ITG *nener,double *trab,
	      ITG *inotr,ITG *ntrans,double *fmpc,char *cbody,ITG *ibody,
	      double *xbody,ITG *nbody,double *xbodyold,ITG *istep,
	      ITG *isolver,ITG *jq,char *output,ITG *mcs,ITG *nkon,
	      ITG *mpcend,ITG *ics,double *cs,ITG *ntie,char *tieset,
	      ITG *idrct,ITG *jmax,double *tmin,double *tmax,
	      double *ctrl,ITG *itpamp,double *tietol,ITG *iit,
	      ITG *ncont,ITG *ne0,double *reltime,double *dtime,
	      double *bcontini,double *bj,double *aux,ITG *iaux,
	      double *bcont,ITG *nev,double *v,
              ITG *nkon0,double *deltmx,double *dtheta,double *theta,
              ITG *iprescribedboundary,ITG *mpcfree,ITG *memmpc_,
              ITG *itietri,ITG *koncont,double *cg,double *straight,
              ITG *iinc,double *vini,
              double *aa,double *bb,double *aanew,double *d,double *z,
	      double *zeta,double *b,double *time0,double *time1,
	      ITG *ipobody,
              double *xforcact,double *xloadact,double *t1act,
              double *xbounact,double *xbodyact,double *cd,double *cv,
              double *ampli,double *dthetaref,double *bjp,double *bp,
              double *cstr,ITG *imddof,
              ITG *nmddof,ITG **ikactcontp,ITG *nactcont,ITG *nactcont_,
              double *aamech,double *bprev,ITG *iprev,ITG *inonlinmpc,
              ITG **ikactmechp,ITG *nactmech,ITG *imdnode,ITG *nmdnode,
              ITG *imdboun,ITG *nmdboun,ITG *imdmpc,ITG *nmdmpc,
              ITG *itp,ITG *inext,
              ITG *imastop,ITG *nslavnode,ITG *islavnode,
              ITG *islavsurf,
              ITG *itiefac,double *areaslav,ITG *iponoels,ITG *inoels,
              double *springarea,ITG *izdof,ITG *nzdof,double *fn,
	      ITG *imastnode,ITG *nmastnode,double *xmastnor,
              double *xstateini,ITG *nslavs,
              ITG *cyclicsymmetry,double *xnoels,ITG *ielas,ITG *ielprop,
              double *prop);
 
void dynboun(double *amta,ITG *namta,ITG *nam,double *ampli,double *time,
             double *ttime,double *dtime,double *xbounold,double *xboun,
             double *xbounact,ITG *iamboun,ITG *nboun,ITG *nodeboun,
             ITG *ndirboun,double *ad,double *au,double *adb,
             double *aub,ITG *icol,ITG *irow,ITG *neq,ITG *nzs,
             double *sigma,double *b,ITG *isolver,
             double *alpham,double *betam,ITG *nzl,
             ITG *init,double *bact,double *bmin,ITG *jq,char *amname,
             double *bv,double *bprev,double *bdiff,
             ITG *nactmech,ITG *icorrect,ITG *iprev);

void FORTRAN(dynresults,(ITG *nk,double *v,ITG *ithermal,ITG *nactdof,
             double *vold,ITG *nodeboun,ITG *ndirboun,double *xboun,
             ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
	     char *labmpc,ITG *nmpc,double *b,double *bp,double *veold,
	     double *dtime,ITG *mi,ITG *imdnode,ITG *nmdnode,ITG *imdboun,
	     ITG *nmdboun,ITG *imdmpc,ITG *nmdmpc,ITG *nmethod,double *time));

void FORTRAN(effectivemodalmass,(ITG *neq,ITG *nactdof,ITG *mi,double *adb,
			double *aub,ITG *jq,ITG *irow,ITG *nev,double *z,
			double *co,ITG *nk));

void electromagnetics(double **co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
	     ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	     ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,
             ITG *nmpc,ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	     ITG **nelemloadp,char **sideloadp,double *xload,
	     ITG *nload,ITG *nactdof,ITG **icolp,ITG *jq,ITG **irowp,
             ITG *neq,ITG *nzl,ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,
             ITG *ikboun,ITG *ilboun,double *elcon,ITG *nelcon,
             double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,
             double *alzero,ITG **ielmatp,ITG **ielorienp,ITG *norien,
             double *orab,ITG *ntmat_,double *t0,double *t1,double *t1old,
	     ITG *ithermal,double *prestr,ITG *iprestr,double **vold,
             ITG *iperturb,double *sti,ITG *nzs,ITG *kode,char *filab,
             ITG *idrct,ITG *jmax,ITG *jout,double *timepar,double *eme,
             double *xbounold,double *xforcold,double *xloadold,
             double *veold,double *accold,char *amname,double *amta,
             ITG *namta,ITG *nam,ITG *iamforc,ITG **iamloadp,ITG *iamt1,
             double *alpha,ITG *iexpl,ITG *iamboun,double *plicon,
             ITG *nplicon,double *plkcon,ITG *nplkcon,
             double **xstatep,ITG *npmat_,ITG *istep,double *ttime,
	     char *matname,double *qaold,ITG *mi,
             ITG *isolver,ITG *ncmat_,ITG *nstate_,ITG *iumat,
             double *cs,ITG *mcs,ITG *nkon,double **ener,ITG *mpcinfo,
             char *output,double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,
             double *physcon,ITG *nflow,double *ctrl,
             char **setp,ITG *nset,ITG **istartsetp,
             ITG **iendsetp,ITG **ialsetp,ITG *nprint,char *prlab,
             char *prset,ITG *nener,ITG *ikforc,ITG *ilforc,double *trab,
             ITG *inotr,ITG *ntrans,double **fmpcp,char *cbody,
             ITG *ibody,double *xbody,ITG *nbody,double *xbodyold,
             ITG *ielprop,double *prop,ITG *ntie,char **tiesetp,
	     ITG *itpamp,ITG *iviewfile,char *jobnamec,double **tietolp,
	     ITG *nslavs,double *thicke,ITG *ics,ITG *nalset,ITG *nmpc_,
	     ITG *nmat,char *typeboun,ITG *iaxial,ITG *nload_);

void FORTRAN(elementpernode,(ITG *iponoel,ITG *inoel,char *lakon,ITG *ipkon,
              ITG *kon,ITG *ne));

void FORTRAN(envtemp,(ITG *itg,ITG *ieg,ITG *ntg,ITG *ntr,char *sideload,
                      ITG *nelemload,ITG *ipkon,ITG *kon,char *lakon,
                      ITG *ielmat,ITG *ne,ITG *nload,
                      ITG *kontri,ITG *ntri,ITG *nloadtr,
                      ITG *nflow,ITG *ndirboun,ITG *nactdog,
                      ITG *nodeboun,ITG *nacteq,
                      ITG *nboun,ITG *ielprop,double *prop,ITG *nteq,
                      double *v,ITG *network,double *physcon,
		      double *shcon,ITG *ntmat_,double *co,
                      double *vold,char *set,ITG *nshcon,
		      double *rhcon,ITG *nrhcon,ITG *mi,ITG *nmpc,
                      ITG *nodempc,ITG *ipompc,char *labmpc,ITG *ikboun,
                      ITG *nasym,ITG *iaxial));

void FORTRAN(equationcheck,(double *ac,ITG *nteq,ITG *nactdog,
                            ITG *itg,ITG *ntg,ITG *nacteq,ITG *network));

void FORTRAN(errorestimator,(double *yi,double *yn,ITG *ipkon,
             ITG *kon,char *lakon,ITG *nk,ITG *ne,ITG *mi,ITG *ielmat,
	     ITG *nterms));

void expand(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	     ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	     ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
             ITG *nmpc,ITG *nodeforc,ITG *ndirforc,double *xforc,
             ITG *nforc,ITG *nelemload,char *sideload,double *xload,
             ITG *nload,ITG *nactdof,ITG *neq,
	     ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,ITG *ilboun,
	     double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	     double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	     ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	     double *t0,ITG *ithermal,double *prestr,ITG *iprestr,
	     double *vold,ITG *iperturb,double *sti,ITG *nzs, 
	     double *adb,double *aub,char *filab,double *eme,
             double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
             double *xstate,ITG *npmat_,char *matname,ITG *mi,
	     ITG *ics,double *cs,ITG *mpcend,ITG *ncmat_,
             ITG *nstate_,ITG *mcs,ITG *nkon,double *ener,
             char *jobnamec,char *output,char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,ITG *nener,double *trab,
             ITG *inotr,ITG *ntrans,double *ttime,double *fmpc,
	     ITG *nev,double **z,ITG *iamboun,double *xbounold,
             ITG *nsectors,ITG *nm,ITG *icol,ITG *irow,ITG *nzl,ITG *nam,
             ITG *ipompcold,ITG *nodempcold,double *coefmpcold,
             char *labmpcold,ITG *nmpcold,double *xloadold,ITG *iamload,
             double *t1old,double *t1,ITG *iamt1,double *xstiff,ITG **icolep,
	     ITG **jqep,ITG **irowep,ITG *isolver,
	     ITG *nzse,double **adbep,double **aubep,ITG *iexpl,ITG *ibody,
	     double *xbody,ITG *nbody,double *cocon,ITG *ncocon,
	     char* tieset,ITG* ntie,ITG *imddof,ITG *nmddof,
	     ITG *imdnode,ITG *nmdnode,ITG *imdboun,ITG *nmdboun,
             ITG *imdmpc,ITG *nmdmpc,ITG **izdofp,ITG *nzdof,ITG *nherm,
	     double *xmr,double *xmi,char *typeboun,ITG *ielprop,double *prop);

void FORTRAN(extrapolate,(double *sti,double *stn,ITG *ipkon,ITG *inum,
             ITG *kon,char *lakon,ITG *nfield,ITG *nk,ITG *ne,ITG *mi,
             ITG *ndim,double *orab,ITG *ielorien,double *co,ITG *iorienglob,
	     char *cflag,ITG *nelemload,ITG *nload,ITG *nodeboun,ITG *nboun,
             ITG *ndirboun,double *vold,ITG *ithermal,ITG *force,
	     ITG *cfd,ITG *ielmat,double *thicke,char *filab));

void FORTRAN(extrapolate_ad_h,(ITG *nface,ITG *ielfa,double *xrlfa,
			       double *ad,double *adfa,
                               double *hel,double *hfa));

void FORTRAN(extrapolate_ad_h_comp,(ITG *nface,ITG *ielfa,double *xrlfa,
			       double *ad,double *adfa,
                               double *hel,double *hfa));

void FORTRAN(extrapolatefluid,(ITG *nk,ITG *iponofa,ITG *inofa,ITG *inum,
             double *vfa,double *v,ITG *ielfa,ITG *ithermal,
             ITG *imach,ITG *ikappa,double *xmach,double *xkappa,
             double *shcon,ITG *nshcon,ITG *ntmat_,ITG *ielmatf,
             double *physcon,ITG *mi));

void FORTRAN(extrapolate_gradtel,(ITG *nface,ITG *ielfa,double *xrlfa,
                    double *gradvel,double *gradvfa));

void FORTRAN(extrapolate_gradvel,(ITG *nface,ITG *ielfa,double *xrlfa,
			       double *gradv,double *gradvfa));

void FORTRAN(extrapolate_pel,(ITG *nface,ITG *ielfa,double *xrlfa,
			      double *vel,double *vfa,ITG *ifabou,
                              double *xboun,ITG *nef));

void FORTRAN(extrapolate_tel,(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
			      double *vfa,ITG *ifabou,double *xboun,ITG *ipnei,
                              ITG *nef));

void FORTRAN(extrapolate_vel,(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
			      double *vfa,ITG *ifabou,double *xboun,ITG *ipnei,
                              ITG *nef));

void FORTRAN(fcrit,(double *time,double *tend,double *aai,double *bbi,
		      double *zetaj,double *dj,double *ddj,
		      double *h1,double *h2,double *h3,double *h4,
                      double *func,double *funcp));

void FORTRAN(findsurface,(ITG *ipoface,ITG *nodface,ITG *ne,ITG *ipkon,ITG *kon,
                     char *lakon,ITG *ntie,char *tieset));

void FORTRAN (flowoutput,(ITG *itg,ITG *ieg,ITG *ntg,ITG *nteq,
			  double *bc,char *lakon,
			  ITG *ntmat_,double *v,double *shcon,ITG *nshcon,
			  ITG *ipkon,ITG *kon,double *co,ITG *nflow,
			  double *dtime,double *ttime,double *time,
			  ITG *ielmat,double *prop,
			  ITG *ielprop,ITG *nactdog,ITG *nacteq,ITG *iin,
			  double *physcon,double *camt,double *camf,double *camp,
			  double *uamt,double *uamf,double *uamp,
			  double *rhcon,ITG *nrhcon,
			  double *vold,char *jobnamef,char *set,ITG *istartset,
                          ITG *iendset,ITG *ialset,ITG *nset,ITG *mi,
                          ITG *iaxial));

void FORTRAN(flowresult,(ITG *ntg,ITG *itg,double *cam,double *vold,
              double *v,
              ITG *nload,char *sideload,ITG *nelemload,
	      double *xloadact,ITG *nactdog,ITG *network,ITG *mi,
	      ITG *ne,ITG *ipkon,char *lakon,ITG *kon));

void FORTRAN(forcesolve,(double *zc,ITG *nev,double *aa,double *bb,
             double *xx,double *eiga,double *eigb,double *eigxx,
	     ITG *iter,double *d,ITG *neq,double *z,ITG *istartnmd,
	     ITG *iendnmd,ITG *nmd,ITG *cyclicsymmetry,ITG *neqact,
	     ITG *igeneralizedforce));

void frd(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne0,
	 double *v,double *stn,ITG *inum,ITG *nmethod,ITG *kode,
	 char *filab,double *een,double *t1,double *fn,double *time,
	 double *epn,ITG *ielmat,char *matname,double *enern,
	 double *xstaten,ITG *nstate_,ITG *istep,ITG *iinc,
	 ITG *ithermal,double *qfn,ITG *mode,ITG *noddiam,
	 double *trab,ITG *inotr,ITG *ntrans,double *orab,
	 ITG *ielorien,ITG *norien,char *description,ITG *ipneigh,
	 ITG *neigh,ITG *mi,double *stx,double *vr,double *vi,
	 double *stnr,double *stni,double *vmax,double *stnmax,
	 ITG *ngraph,double *veold,double *ener,ITG *ne,double *cs,
	 char *set,ITG *nset,ITG *istartset,ITG *iendset,ITG *ialset,
	 double *eenmax,double *fnr,double *fni,double *emn,
	 double *thicke,char *jobnamec,char *output,double *qfx,
         double *cdn,ITG *mortar,double *cdnr,double *cdni,ITG *nmat);

void frdcyc(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,double *v,
	    double *stn,ITG *inum,ITG *nmethod,ITG *kode,char *filab,
	    double *een,double *t1,double *fn,double *time,double *epn,
	    ITG *ielmat,char *matname,double *cs,ITG *mcs,ITG *nkon,
	    double *enern,double *xstaten,ITG *nstate_,ITG *istep,
            ITG *iinc,ITG *iperturb,double *ener,ITG *mi,char *output,
            ITG *ithermal,double *qfn,ITG *ialset,ITG *istartset,
            ITG *iendset,double *trab,ITG *inotr,ITG *ntrans,double *orab,
	    ITG *ielorien,ITG *norien,double *sti,double *veold,ITG *noddiam,
            char *set,ITG *nset,double *emn,double *thicke,char *jobnamec,
            ITG *ne0,double *cdn,ITG *mortar,ITG *nmat);

void frd_se(double *co,ITG *nk,double *stn,ITG *inum,ITG *nmethod,
         ITG *kode,char *filab,double *fn,double *time,ITG *nstate_,
	 ITG *istep,ITG *iinc,ITG *mode,ITG *noddiam,char *description,
	 ITG *mi,ITG *ngraph,ITG *ne,double *cs,char *set,ITG *nset,
	 ITG *istartset,ITG *iendset,ITG *ialset,double *thicke,
	 char *jobnamec,char *output,double *dgdxtotglob,ITG *numobject);

void FORTRAN(frdfluid,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
             ITG *ne,double *v,ITG *kode,double *time,
             ITG *ielmatf,char *matname,char *filab,
	     ITG *inum,ITG *ntrans,ITG *inotr,double *trab,ITG *mi,
	     ITG *istep,double *stn,double *qfn,ITG *nactdofinv,
	     double *xmach,double *xkappa,double *physcon));

void frdheader(ITG *icounter,double *oner,double *time,double *pi,
	       ITG *noddiam,double *cs,ITG *null,ITG *mode,
	       ITG *noutloc,char *description,ITG *kode,ITG *nmethod,
               FILE *f1,char *output,ITG *istep,ITG *iinc);

void FORTRAN(frditeration,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
             ITG *ne,double *v,double *time,ITG *ielmat,char *matname,
	     ITG *mi,ITG *istep,ITG *iinc,ITG *ithermal));

void frdselect(double *field1,double *field2,ITG *iset,ITG *nkcoords,ITG *inum,
     char *m1,ITG *istartset,ITG *iendset,ITG *ialset,ITG *ngraph,ITG *ncomp,
     ITG *ifield,ITG *icomp,ITG *nfield,ITG *iselect,char *m2,FILE *f1,
     char *output,char *m3);

void frdset(char *filabl,char *set,ITG *iset,ITG *istartset,ITG *iendset,
	    ITG *ialset,ITG *inum,ITG *noutloc,ITG *nout,ITG *nset,
	    ITG *noutmin,ITG *noutplus,ITG *iselect,ITG *ngraph);

void frdvector(double *v,ITG *iset,ITG *ntrans,char * filabl,ITG *nkcoords,
               ITG *inum,char *m1,ITG *inotr,double *trab,double *co,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *mi,ITG *ngraph,
               FILE *f1,char *output,char *m3);

void FORTRAN(frictionheating,(ITG *ne0,ITG *ne,ITG *ipkon,char *lakon,ITG *ielmat,
                     ITG *mi,double *elcon,ITG *ncmat_,ITG *ntmat_,
                     ITG *kon,ITG *islavsurf,double *pmastsurf,
                     double *springarea,double *co,double *vold,
                     double *veold,double *pslavsurf,double *xload,
                     ITG *nload,ITG *nload_,ITG *nelemload,ITG *iamload,
		     ITG *idefload,char *sideload,double *stx,ITG *nam));

void FORTRAN(fsub,(double *time,double *tend,double *aai,double *bbi,
		   double *ddj,double *h1,double *h2,double *h3,double *h4,
                   double *func,double *funcp));

void FORTRAN(fsuper,(double *time,double *tend,double *aai,double *bbi,
		       double *h1,double *h2,double *h3,double *h4,
		       double *h5,double *h6,double *func,double *funcp));

void FORTRAN(gasmechbc,(double *vold,ITG *nload,char *sideload,
			ITG *nelemload,double *xload,ITG *mi));

void FORTRAN(genadvecelem,(ITG *inodesd,ITG *ipkon,ITG *ne,char *lakon,
             ITG *kon,ITG *nload,char *sideload,ITG *nelemload,ITG *nkon));

void FORTRAN(gencontelem_f2f,(char *tieset,ITG *ntie,ITG *itietri,ITG *ne,
             ITG *ipkon,ITG *kon,char *lakon,double *cg,double *straight,
             ITG *ifree,ITG *koncont,double *co,double *vold,double *xo,
             double *yo,double *zo,double *x,double *y,double *z,ITG *nx,
             ITG *ny,ITG *nz,ITG *ielmat,double *elcon,ITG *istep,ITG *iinc,
             ITG *iit,ITG *ncmat_,ITG *ntmat_,ITG *mi,ITG *imastop,
             ITG *islavsurf,ITG *itiefac,double *springarea,double *tietol,
             double *reltime,char *filab,ITG *nasym,
	     double *pslavsurf,double *pmastsurf,double *clearini,
	     double *theta,double *xstateini,double *xstate,ITG *nstate_,
	     ITG *ne0,ITG *icutb,ITG *ialeatoric,ITG *nmethod));

void FORTRAN(gencontelem_n2f,(char *tieset,ITG *ntie,ITG *itietri,ITG *ne,
     ITG *ipkon,ITG *kon,char *lakon,
     double *cg,double *straight,ITG *ifree,ITG *koncont,
     double *co,double *vold,double *xo,double *yo,double *zo,
     double *x,double *y,double *z,ITG *nx,ITG *ny,ITG *nz,
     ITG *ielmat,double *elcon,ITG *istep,ITG *iinc,ITG *iit,
     ITG *ncmat_,ITG *ntmat_,
     ITG *nmethod,ITG *mi,ITG *imastop,ITG *nslavnode,
     ITG *islavnode,ITG *islavsurf,ITG *itiefac,double *areaslav,
     ITG *iponoels,ITG *inoels,double *springarea,
     char *set,ITG *nset,ITG *istartset,ITG *iendset,ITG *ialset,
     double *tietol,double *reltime,
     char* filab,ITG *nasym,double *xnoels,ITG *icutb));

void FORTRAN(generateeminterfaces,(ITG *istartset,ITG *iendset,
	     ITG *ialset,ITG *iactive,ITG *ipkon,char *lakon,ITG *kon,
	     ITG *ikmpc,ITG *nmpc,ITG *nafaces));

void FORTRAN(generatetet,(ITG *kontet,ITG *ifatet,ITG *netet,
             ITG *inodfa,ITG *ifreefa,double *planfa,ITG *ipofa,
             ITG *nodes,double *cotet));

void  FORTRAN(gennactdofinv,(ITG *nactdof,ITG *nactdofinv,ITG *nk,
       ITG *mi,ITG *nodorig,ITG *ipkon,char *lakon,ITG *kon,ITG *ne));

void FORTRAN(gentiedmpc,(char *tieset,ITG *ntie,ITG *itietri,
          ITG *ipkon,ITG *kon,char *lakon,char *set,ITG *istartset,
          ITG *iendset,ITG *ialset,double *cg,double *straight,
	  ITG *koncont,double *co,double *xo,double *yo,double *zo,
          double *x,double *y,double *z,ITG *nx,ITG *ny,ITG *nz,ITG *nset,
          ITG *ifaceslave,ITG *istartfield,ITG *iendfield,ITG *ifield,
          ITG *ipompc,ITG *nodempc,double *coefmpc,ITG *nmpc,ITG *nmpc_,
          ITG *mpcfree,ITG *ikmpc,ITG *ilmpc,char *labmpc,ITG *ithermal,
	  double *tietol,ITG *icfd,ITG *ncont,ITG *imastop,ITG *ikboun,
	  ITG *nboun,char *kind));

void FORTRAN(geomview,(double *vold,double *co,double *pmid,double *e1,
             double *e2,double *e3,ITG *kontri,double *area,double *cs,
             ITG *mcs,ITG *inocs,ITG *ntrit,ITG *nk,ITG *mi,double *sidemean));

void FORTRAN(getdesiinfo,(char *set,ITG *istartset,ITG *iendset,ITG *ialset,
            ITG *nset,ITG *mi,ITG *nactdof,ITG *ndesi,
            ITG *ndirdesi,ITG *nodedesi,ITG *ntie,char *tieset));  

void getglobalresults (char *jobnamec,ITG **integerglobp,double **doubleglobp,
                       ITG *nboun,ITG *iamboun,double *xboun,ITG *nload,
                       char *sideload,ITG *iamload,ITG *iglob,ITG *nforc,
                       ITG *iamforc,double *xforc,ITG *ithermal,ITG *nk,
                       double *t1,ITG *iamt1);

ITG getSystemCPUs();;

void FORTRAN(identamta,(double *amta,double *reftime,ITG *istart,ITG *iend,
               ITG *id));

void FORTRAN(identifytiedface,(char *tieset,ITG *ntie,char *set,ITG *nset,
			       ITG *faceslave,char *kind));

void FORTRAN(includefilename,(char *buff,char *includefn,ITG *lincludefn));

void FORTRAN(inicalcbody,(ITG *nef,double *body,ITG *ipobody,ITG *ibody,
		          double *xbody,double *coel,double *vel,char *lakon,
			  ITG *nactdohinv,ITG *icent));

void inicont(ITG* nk,ITG *ncont,ITG *ntie,char *tieset,ITG *nset,char *set,
             ITG *istartset,ITG *iendset,ITG *ialset,ITG **itietrip,
	     char *lakon,ITG *ipkon,ITG *kon,ITG **koncontp,
             ITG *ncone,double *tietol,ITG *ismallsliding,ITG **itiefacp,
	     ITG **islavsurfp,ITG **islavnodep,ITG **imastnodep,
	     ITG **nslavnodep,ITG **nmastnodep,ITG *mortar,
	     ITG **imastopp,ITG *nkon,ITG **iponoels,ITG **inoelsp,
             ITG **ipep,ITG **imep,ITG *ne,ITG *ifacecount,
             ITG *nmpc,ITG *mpcfree,ITG *memmpc_,
             ITG **ipompcp,char **labmpcp,ITG **ikmpcp,ITG **ilmpcp,
             double **fmpcp,ITG **nodempcp,double **coefmpcp,
             ITG *iperturb,ITG *ikboun,ITG *nboun,double *co,ITG *istep,
             double **xnoelsp);

void FORTRAN(init,(ITG *nktet,ITG *inodfa,ITG *ipofa,ITG *netet_));

void FORTRAN(initialcfd,(ITG *ne,ITG *ipkon,ITG *kon,char *lakon,
             double *co,double *coel,double *cofa,ITG *nface,ITG *ielfa,
             double *area,ITG *ipnei,ITG *neiel,double *xxn,double *xxi,
	     double *xle,double *xlen,double *xlet,double *xrlfa,double *cosa,
	     double *volume,ITG *neifa,double *xxj,double *cosb,double *vel,
             double *dmin));

void FORTRAN(initialchannel,(ITG *itg,ITG *ieg,ITG *ntg,double *ac,double *bc,
                         char *lakon,double *v,ITG * ipkon,ITG *kon,
                         ITG *nflow,ITG *ikboun,ITG *nboun,double *prop,
                         ITG *ielprop,ITG *nactdog,ITG *ndirboun,
                         ITG *nodeboun,double *xbounact,ITG *ielmat,
                         ITG *ntmat_,double *shcon,ITG *nshcon,
                         double *physcon,ITG *ipiv,ITG *nteq,
                         double *rhcon,ITG *nrhcon,ITG *ipobody,ITG *ibody,
                         double *xbody,double *co,ITG *nbody,ITG *network,
                         ITG *iin_abs,double *vold,char *set,ITG *istep,
                         ITG *iit,ITG *mi,ITG *ineighe,ITG *ilboun));

void FORTRAN(initialnet,(ITG *itg,ITG *ieg,ITG *ntg,double *ac,double *bc,
                         char *lakon,double *v,ITG * ipkon,ITG *kon,
                         ITG *nflow,ITG *ikboun,ITG *nboun,double *prop,
                         ITG *ielprop,ITG *nactdog,ITG *ndirboun,
                         ITG *nodeboun,double *xbounact,ITG *ielmat,
                         ITG *ntmat_,double *shcon,ITG *nshcon,
                         double *physcon,ITG *ipiv,ITG *nteq,
                         double *rhcon,ITG *nrhcon,ITG *ipobody,ITG *ibody,
                         double *xbody,double *co,ITG *nbody,ITG *network,
                         ITG *iin_abs,double *vold,char *set,ITG *istep,
                         ITG *iit,ITG *mi,ITG *ineighe,ITG *ilboun,
                         ITG *channel,ITG *iaxial,ITG *nmpc,char *labmpc,
                         ITG *ipompc,ITG *nodempc,double *coefmpc));

void insert(ITG *ipointer,ITG **mast1p,ITG **mast2p,ITG *i1,
	    ITG *i2,ITG *ifree,ITG *nzs_);

void insertrad(ITG *ipointer,ITG **mast1p,ITG **mast2p,ITG *i1,
	    ITG *i2,ITG *ifree,ITG *nzs_);

void FORTRAN(integral_boundary,(double *sumfix,double *sumfree,ITG *ifaext,
				ITG *nfaext,ITG *ielfa,ITG *ifabou,double *vfa,ITG *ipnei,
                                double *xxn));

void FORTRAN(interpolatestate,(ITG *ne,ITG *ipkon,ITG *kon,char *lakon,
             ITG *ne0,ITG *mi,double *xstate,
             double *pslavsurf,ITG *nstate_,double *xstateini,
             ITG *islavsurf,ITG *islavsurfold,
	     double *pslavsurfold,char *tieset,ITG *ntie,ITG *itiefac));

void FORTRAN(islavactive,(char *tieset,ITG *ntie,ITG *itietri,double *cg,
              double *straight,double *co,double *vold,double *xo,
              double *yo,double *zo,double *x,double *y,double *z,
              ITG *nx,ITG *ny,ITG *nz,ITG *mi,ITG *imastop,ITG *nslavnode,
              ITG *islavnode,ITG *islavact));

void FORTRAN(isortid,(ITG *ix,double *dy,ITG *n,ITG *kflag));

void FORTRAN(isortii,(ITG *ix,ITG *iy,ITG *n,ITG *kflag));

void FORTRAN(isortiid,(ITG *ix,ITG *iy,double *dy,ITG *n,ITG *kflag));

void FORTRAN(isortiddc,(ITG *ix,double *dy1,double *dy2,char *cy,ITG *n,
                         ITG *kflag));

void FORTRAN(isortiiddc,(ITG *ix1,ITG *ix2,double *dy1,double *dy2,
                         char *cy,ITG *n,ITG *kflag));

void FORTRAN(jouleheating,(ITG *ipkon,char *lakon,ITG *kon,double *co,
             double *elcon,ITG *nelcon,ITG *mi,ITG *ne,double *sti,
             ITG *ielmat,ITG *nelemload,char *sideload,double *xload,
             ITG *nload,ITG *nload_,ITG *iamload,ITG *nam,ITG *idefload,
             ITG *ncmat_,ITG *ntmat_,double *alcon,ITG *nalcon,ITG *ithermal,
             double *vold,double *t1));

void FORTRAN(keystart,(ITG *ifreeinp,ITG *ipoinp,ITG *inp,char *name,
           ITG *iline,ITG *ikey));
  
void linstatic(double *co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
	     ITG *ne,
	     ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	     ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
             ITG *nmpc,
	     ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	     ITG *nelemload,char *sideload,double *xload,
	     ITG *nload,ITG *nactdof,
	     ITG **icolp,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
	     ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	     ITG *ilboun,
	     double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	     double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
	     ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	     double *t0,double *t1,double *t1old,
	     ITG *ithermal,double *prestr,ITG *iprestr,
	     double *vold,ITG *iperturb,double *sti,ITG *nzs,
	     ITG *kode,char *filab,double *eme,
             ITG *iexpl,double *plicon,ITG *nplicon,double *plkcon,
             ITG *nplkcon,
             double **xstatep,ITG *npmat_,char *matname,ITG *isolver,
	     ITG *mi,ITG *ncmat_,ITG *nstate_,double *cs,ITG *mcs,
             ITG *nkon,double **enerp,double *xbounold,
	     double *xforcold,double *xloadold,
             char *amname,double *amta,ITG *namta,
             ITG *nam,ITG *iamforc,ITG *iamload,
             ITG *iamt1,ITG *iamboun,double *ttime,char *output,
             char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,ITG *nener,double *trab,
             ITG *inotr,ITG *ntrans,double *fmpc,char *cbody,ITG *ibody,
	     double *xbody,ITG *nbody,double *xbodyold,double *timepar,
	     double *thicke,char *jobnamec,char *tieset,ITG *ntie,
	     ITG *istep,ITG *nmat,ITG *ielprop,double *prop,char *typeboun,
	     ITG *mortar,ITG *mpcinfo,double *tietol,ITG *ics,ITG *icontact);

void FORTRAN(mafillcorio,(double *co,ITG *nk,ITG *kon,ITG *ipkon,
               char *lakon,
	       ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
	       ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
	       ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
	       ITG *nbody,double *cgr,
	       double *ad,double *au,ITG *nactdof,
	       ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,
	       double *prestr,ITG *iprestr,double *vold,
	       ITG *iperturb,double *sti,ITG *nzs,double *stx,
	       double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
	       ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,double *ttime,double *time,
	       ITG *istep,ITG *kinc,ITG *ibody,ITG *ielprop,double *prop));

void FORTRAN(mafilldm,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
	       ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
	       ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
	       ITG *nbody,double *cgr,
	       double *ad,double *au,ITG *nactdof,
	       ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,
	       double *prestr,ITG *iprestr,double *vold,
	       ITG *iperturb,double *sti,ITG *nzs,double *stx,
	       double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
	       ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *ibody,double *clearini,
               ITG *mortar,double *springarea,double *pslavsurf,
               double *pmastsurf,double *reltime,ITG *nasym));

void FORTRAN(mafillem,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
	       ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
	       ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
	       ITG *nbody,double *cgr,
	       double *ad,double *au,double *bb,ITG *nactdof,
	       ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,
	       double *prestr,ITG *iprestr,double *vold,
	       ITG *iperturb,double *sti,ITG *nzs,double *stx,
	       double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
	       ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
	       double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
	       double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
	       ITG *nasym,ITG *iactive,double *h0,double *pslavsurf,
	       double *pmastsurf,ITG *mortar,double *clearini,
	       ITG *ielprop,double *prop));

void FORTRAN(mafillnet,(ITG *itg,ITG *ieg,ITG *ntg,
			double *ac,ITG *nload,char *sideload,
			ITG *nelemload,double *xloadact,char *lakon,
			ITG *ntmat_,double *v,double *shcon,ITG *nshcon,
			ITG *ipkon,ITG *kon,double *co,ITG *nflow,
			ITG *iinc,ITG *istep,
			double *dtime,double *ttime,double *time,
			ITG *ielmat,ITG *nteq,double *prop,
			ITG *ielprop,ITG *nactdog,ITG *nacteq,
			double *physcon,double *rhcon,ITG *nrhcon,
			ITG *ipobody,ITG *ibody,double *xbody,ITG *nbody,
			double *vold,double *xloadold,double *reltime,
			ITG *nmethod,char *set,ITG *mi,ITG *nmpc,
                        ITG *nodempc,ITG *ipompc,double *coefmpc,
                        char *labmpc,ITG *iaxial));

void FORTRAN(mafillp,(ITG *ne,char *lakonf,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
	     double *xle,double *b,double *xxn,
	     ITG *neq,ITG *nzs,double *hfa,double *gradpel,
	     double *bp,double *xxi,ITG *neij,double *xlen,double *cosb,
             ITG *nefa,ITG *nefb));

void FORTRAN(mafillpbc,(ITG *nef,double *au,double *ad,ITG *jq,ITG *irow,
			double *b,ITG *iatleastonepressurebc,ITG *nzs));

void FORTRAN(mafillpcomp,(ITG *ne,char *lakonf,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
	     double *xle,double *b,double *xxn,ITG *neq,
	     ITG *nzs,double *hfa,double *gradpel,
	     double *bp,double *xxi,ITG *neij,double *xlen,double *cosb,
             ITG *ielmatf,ITG *mi,double *a1,double *a2,double *a3,double *velo,
             double *veloo,double *dtimef,double *shcon,ITG *ntmat_,double *vel,
	     ITG *nactdohinv,double *xrlfa,double *flux,ITG *nefa,ITG *nefb));

void mafillpcompmain(ITG *ne,char *lakonf,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
	     double *xle,double *b,double *xxn,ITG *neq,
	     ITG *nzs,double *hfa,double *gradpel,
	     double *bp,double *xxi,ITG *neij,double *xlen,double *cosb,
             ITG *ielmatf,ITG *mi,double *a1,double *a2,double *a3,double *velo,
             double *veloo,double *dtimef,double *shcon,ITG *ntmat_,double *vel,
	     ITG *nactdohinv,double *xrlfa,double *flux);

void *mafillpcompmt(ITG *i);

void mafillpmain(ITG *ne,char *lakonf,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
	     double *xle,double *b,double *xxn,
	     ITG *neq,ITG *nzs,double *hfa,double *gradpel,
	     double *bp,double *xxi,ITG *neij,double *xlen,double *cosb,
	     ITG *iatleastonepressurebc);

void *mafillpmt(ITG *i);

void FORTRAN(mafillsm,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
	       ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
	       ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
	       ITG *nbody,double *cgr,
	       double *ad,double *au,double *bb,ITG *nactdof,
	       ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,
	       double *prestr,ITG *iprestr,double *vold,
	       ITG *iperturb,double *sti,ITG *nzs,double *stx,
	       double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
	       ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
	       double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
	       double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
	       ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
	       double *clearini,ITG *ielprop,double *prop,ITG *ne0,
	       double *fnext,ITG *nea,ITG *neb,ITG *kscale));

void FORTRAN(mafillsmse,(double *co,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,ITG *ipompc,ITG *nodempc,double *coefmpc,
	       ITG *nmpc,ITG *nelemload,char *sideload,double *xload,
	       ITG *nload,double *xbody,ITG *ipobody,ITG *nbody,
	       double *cgr,ITG *nactdof,ITG *neq,ITG *nmethod,ITG *ikmpc,
	       ITG *ilmpc,double *elcon,ITG *nelcon,double *rhcon,
	       ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
	       ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,
	       ITG *ntmat_,double *t0,double *t1,ITG *ithermal,
	       ITG *iprestr,double *vold,ITG *iperturb,double *sti,
	       double *stx,ITG *iexpl,double *plicon,ITG *nplicon,
               double *plkcon,ITG *nplkcon,double *xstiff,ITG *npmat_,
	       double *dtime,char *matname,ITG *mi,ITG *ncmat_,ITG *mass,
               ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *ttime,double *time,
               ITG *istep,ITG *iinc,ITG *coriolis,ITG *ibody,
	       double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
	       double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
	       ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
	       double *clearini,ITG *ielprop,double *prop,ITG *ne0,
	       ITG *nea,ITG *neb,double *distmin,ITG *ndesi,
	       ITG *nodedesi,ITG *ndirdesi,double *dfextminds));

void *mafillsmmt(ITG *i);

void *mafillsmsemt(ITG *i);

void mafillsmmain(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
	       ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
	       ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
	       ITG *nbody,double *cgr,
	       double *ad,double *au,double *bb,ITG *nactdof,
	       ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,
	       double *prestr,ITG *iprestr,double *vold,
	       ITG *iperturb,double *sti,ITG *nzs,double *stx,
	       double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
	       ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
	       double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
	       double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
	       ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
	       double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               double *fnext,ITG *kscale);

void mafillsmmain_se(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
	       ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
	       ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
	       ITG *nbody,double *cgr,
	       double *ad,double *au,ITG *nactdof,
	       ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,
	       double *prestr,ITG *iprestr,double *vold,
	       ITG *iperturb,double *sti,ITG *nzs,double *stx,
	       double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
	       ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
	       double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
	       double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
	       ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
	       double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               double *fnext,double *distmin,ITG *ndesi,ITG *nodedesi,
	       ITG *ndirdesi,double *dfextminds);

void FORTRAN(mafillsmas,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
	       ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
	       ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
	       ITG *nbody,double *cgr,
	       double *ad,double *au,double *bb,ITG *nactdof,
	       ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,
	       double *prestr,ITG *iprestr,double *vold,
	       ITG *iperturb,double *sti,ITG *nzs,double *stx,
	       double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
	       ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
	       double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
	       double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
	       ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               ITG *kscale));

void FORTRAN(mafillsmas1,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
	       ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
	       ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
	       ITG *nbody,double *cgr,
	       double *ad,double *au,double *bb,ITG *nactdof,
	       ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,
	       double *prestr,ITG *iprestr,double *vold,
	       ITG *iperturb,double *sti,ITG *nzs,double *stx,
	       double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
	       ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
	       double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
	       double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
	       ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
	       double *clearini,ITG *ielprop,double *prop,ITG *ne0,
	       ITG *nea,ITG *neb,ITG *kscale));

void mafillsmasmain(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
	       ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
	       ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
	       ITG *nbody,double *cgr,
	       double *ad,double *au,double *bb,ITG *nactdof,
	       ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,
	       double *prestr,ITG *iprestr,double *vold,
	       ITG *iperturb,double *sti,ITG *nzs,double *stx,
	       double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
	       ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
	       double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
	       double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
	       ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               ITG *kscale);

void *mafillsmasmt(ITG *i);

void FORTRAN(mafillsmcs,(double *co,ITG *nk,ITG *kon,ITG *ipkon,
               char *lakon,
	       ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
	       ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
	       ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
	       ITG *nbody,double *cgr,
	       double *ad,double *au,double *bb,ITG *nactdof,
	       ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,
	       double *prestr,ITG *iprestr,double *vold,
	       ITG *iperturb,double *sti,ITG *nzs,double *stx,
	       double *adb,double *aub,ITG *iexpl,double *plicon,
               ITG *nplicon,double *plkcon,ITG *nplkcon,double *xstiff,
	       ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ics,double *cs,ITG *nm,ITG *ncmat_,char *labmpc,
               ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,ITG *mcs,ITG *coriolis,ITG *ibody,
	       double *xloadold,double *reltime,ITG *ielcs,double *veold,
	       double *springarea,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
	       ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               ITG *kscale));

void FORTRAN(mafillsmcsas,(double *co,ITG *nk,ITG *kon,ITG *ipkon,
               char *lakon,
	       ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
	       ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
	       ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
	       ITG *nbody,double *cgr,
	       double *ad,double *au,double *bb,ITG *nactdof,
	       ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,
	       double *prestr,ITG *iprestr,double *vold,
	       ITG *iperturb,double *sti,ITG *nzs,double *stx,
	       double *adb,double *aub,ITG *iexpl,double *plicon,
               ITG *nplicon,double *plkcon,ITG *nplkcon,double *xstiff,
	       ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ics,double *cs,ITG *nm,ITG *ncmat_,char *labmpc,
               ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,ITG *mcs,ITG *coriolis,ITG *ibody,
	       double *xloadold,double *reltime,ITG *ielcs,double *veold,
	       double *springarea,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
	       ITG *nasym,ITG *nstate_,double *xstateini,double *xstate,
	       double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               ITG *kscale));

void FORTRAN(mafillsmforc,(ITG *nforc,ITG *ndirforc,ITG *nodeforc,
             double *xforc,ITG *nactdof,double *fext,ITG *nmpc,ITG *ipompc,
             ITG *nodempc,ITG *ikmpc,ITG *ilmpc,double *coefmpc,ITG *mi,
	     ITG *rhsi,double *fnext,ITG *nmethod));

void FORTRAN(mafillt,(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
	       double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umel,double *xlet,
               double *xle,double *gradtfa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
	       double *gradvel,double *xload,double *gammat,double *xrlfa,
	       double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
	       double *flux,ITG *nefa,ITG *nefb));

void FORTRAN(mafilltcomp,(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
	       double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umel,double *xlet,
               double *xle,double *gradtfa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
	       double *gradvel,double *xload,double *gammat,double *xrlfa,
	       double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux,ITG *nefa,ITG *nefb));

void mafilltcompmain(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
	       double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umel,double *xlet,
               double *xle,double *gradtfa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
	       double *gradvel,double *xload,double *gammat,double *xrlfa,
	       double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux);

void *mafilltcompmt(ITG *i);

void mafilltmain(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
	       double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umel,double *xlet,
               double *xle,double *gradtfa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
	       double *gradvel,double *xload,double *gammat,double *xrlfa,
	       double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
	       double *flux);

void *mafilltmt(ITG *i);

void FORTRAN(mafillv,(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,
             double *vfa,double *xxn,double *area,double *au,double *ad,
             ITG *jq,ITG *irow,ITG *nzs,double *b,double *vel,double *cosa,
             double *umfa,double *xlet,double *xle,double *gradvfa,
             double *xxi,double *body,double *volume,
	     ITG *ielfa,char *lakonf,ITG *ifabou,ITG *nbody,ITG *neq,
	     double *dtimef,double *velo,double *veloo,
	     double *sel,double *xrlfa,double *gamma,double *xxj,
	     ITG *nactdohinv,double *a1,double *a2,double *a3,double *flux,
	     ITG *nefa,ITG *nefb));

void FORTRAN(mafillvcomp,(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,
             double *vfa,double *xxn,double *area,double *au,double *ad,
             ITG *jq,ITG *irow,ITG *nzs,double *b,double *vel,double *cosa,
             double *umfa,double *xlet,double *xle,double *gradvfa,
	     double *xxi,double *body,double *volume,
	     ITG *ielfa,char *lakonf,ITG *ifabou,ITG *nbody,ITG *neq,
	     double *dtimef,double *velo,double *veloo,
	     double *sel,double *xrlfa,double *gamma,double *xxj,
	     ITG *nactdohinv,double *a1,double *a2,double *a3,double *flux,
	     ITG *nefa,ITG *nefb));

void mafillvcompmain(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,
             double *vfa,double *xxn,double *area,double *au,double *ad,
             ITG *jq,ITG *irow,ITG *nzs,double *b,double *vel,double *cosa,
             double *umfa,double *xlet,double *xle,double *gradvfa,
	     double *xxi,double *body,double *volume,
	     ITG *ielfa,char *lakonf,ITG *ifabou,ITG *nbody,ITG *neq,
	     double *dtimef,double *velo,double *veloo,
	     double *sel,double *xrlfa,double *gamma,double *xxj,
	     ITG *nactdohinv,double *a1,double *a2,double *a3,double *flux);

void *mafillvcompmt(ITG *i);

void mafillvmain(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,
             double *vfa,double *xxn,double *area,double *au,double *ad,
             ITG *jq,ITG *irow,ITG *nzs,double *b,double *vel,double *cosa,
             double *umfa,double *xlet,double *xle,double *gradvfa,
             double *xxi,double *body,double *volume,
	     ITG *ielfa,char *lakonf,ITG *ifabou,ITG *nbody,ITG *neq,
	     double *dtimef,double *velo,double *veloo,
	     double *sel,double *xrlfa,double *gamma,double *xxj,
	     ITG *nactdohinv,double *a1,double *a2,double *a3,double *flux);

void *mafillvmt(ITG *i);

void mastruct(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
	      ITG *nodeboun,ITG *ndirboun,ITG *nboun,ITG *ipompc,
	      ITG *nodempc,ITG *nmpc,ITG *nactdof,ITG *icol,
	      ITG *jq,ITG **mast1p,ITG **irowp,ITG *isolver,ITG *neq,
	      ITG *ikmpc,ITG *ilmpc,ITG *ipointer,ITG *nzs,ITG *nmethod,
              ITG *ithermal,ITG *ikboun,ITG *ilboun,ITG *iperturb,
              ITG *mi,ITG *mortar,char *typeboun,char *labmpc,
              ITG *iit,ITG *icascade);

void mastructcs(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,ITG *nodeboun,
	       ITG *ndirboun,ITG *nboun,ITG *ipompc,ITG *nodempc,
	       ITG *nmpc,ITG *nactdof,ITG *icol,ITG *jq,ITG **mast1p,
	       ITG **irowp,ITG *isolver,ITG *neq,
	       ITG *ikmpc,ITG *ilmpc,ITG *ipointer,
	       ITG *nzs,ITG *nmethod,ITG *ics,double *cs,
	       char *labmpc,ITG *mcs,ITG *mi,ITG *mortar);

void mastructem(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
	      ITG *nodeboun,ITG *ndirboun,ITG *nboun,ITG *ipompc,
	      ITG *nodempc,ITG *nmpc,ITG *nactdof,ITG *icol,
	      ITG *jq,ITG **mast1p,ITG **irowp,ITG *isolver,ITG *neq,
	      ITG *ikmpc,ITG *ilmpc,ITG *ipointer,ITG *nzs,
	      ITG *ithermal,ITG *mi,ITG *ielmat,double *elcon,ITG *ncmat_,
	      ITG *ntmat_,ITG *inomat);

void mastructf(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
	       ITG *icol,ITG *jq,ITG **mast1p,ITG **irowp,
	       ITG *isolver,ITG *neq,ITG *ipointer,ITG *nzs,
               ITG *ipnei,ITG *ineiel,ITG *mi);

void mastructrad(ITG *ntr,ITG *nloadtr,char *sideload,ITG *ipointerrad,
              ITG **mast1radp,ITG **irowradp,ITG *nzsrad,
	      ITG *jqrad,ITG *icolrad);

void FORTRAN(materialdata_cfd,(ITG *nef,double *vel,double *shcon,
	     ITG *nshcon,ITG *ielmat,ITG *ntmat_,ITG *mi,double *cvel,
             double *vfa,double *cocon,ITG *ncocon,double *physcon,
             double *cvfa,ITG *ithermal,ITG *nface,double *umel,
	     double *umfa,ITG *ielfa,double *hcfa,double *rhcon,
             ITG *nrhcon));

void FORTRAN(materialdata_cfd_comp,(ITG *nef,double *vel,double *shcon,
	     ITG *nshcon,ITG *ielmat,ITG *ntmat_,ITG *mi,double *cvel,
             double *vfa,double *cocon,ITG *ncocon,double *physcon,
             double *cvfa,ITG *ithermal,ITG *nface,double *umel,
             double *umfa,ITG *ielfa,double *hcfa));

void FORTRAN(meannode,(ITG *nk,ITG *inum,double *v));

void FORTRAN(mpcrem,(ITG *i,ITG *mpcfree,ITG *nodempc,ITG *nmpc,ITG *ikmpc,
                     ITG *ilmpc,char *labmpc,double *coefmpc,ITG *ipompc));

void FORTRAN(mult,(double *matrix,double *trans,ITG *n));

void FORTRAN(negativepressure,(ITG *ne0,ITG *ne,ITG *mi,double *stx,
			       double *pressureratio));

void FORTRAN(networkinum,(ITG *ipkon,ITG *inum,ITG *kon,char *lakon,
       ITG *ne,ITG *itg,ITG *ntg));

void FORTRAN(nident,(ITG *x,ITG *px,ITG *n,ITG *id));

void FORTRAN(nidentll,(long long *x,long long *px,ITG *n,ITG *id));

void FORTRAN(nodestiedface,(char *tieset,ITG *ntie,ITG *ipkon,ITG *kon,
       char *lakon,char *set,ITG *istartset,ITG *iendset,ITG *ialset,
       ITG *nset,ITG *faceslave,ITG *istartfield,ITG *iendfield,
       ITG *ifield,ITG *nconf,ITG *ncone,char *kind));

void nonlingeo(double **co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
	     ITG *ne,
	     ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	     ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,
             ITG *nmpc,
	     ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	     ITG **nelemloadp,char **sideloadp,double *xload,
	     ITG *nload,ITG *nactdof,
	     ITG **icolp,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
	     ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,ITG *ikboun,
	     ITG *ilboun,
	     double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	     double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
	     ITG **ielorienp,ITG *norien,double *orab,ITG *ntmat_,
	     double *t0,double *t1,double *t1old,
	     ITG *ithermal,double *prestr,ITG *iprestr,
	     double **vold,ITG *iperturb,double *sti,ITG *nzs, 
	     ITG *kode,char *filab,ITG *idrct,
	     ITG *jmax,ITG *jout,double *timepar,
	     double *eme,double *xbounold,
	     double *xforcold,double *xloadold,
             double *veold,double *accold,
             char *amname,double *amta,ITG *namta,ITG *nam,
             ITG *iamforc,ITG **iamloadp,
             ITG *iamt1,double *alpha,ITG *iexpl,
	     ITG *iamboun,double *plicon,ITG *nplicon,double *plkcon,
	     ITG *nplkcon,
             double **xstatep,ITG *npmat_,ITG *istep,double *ttime,
	     char *matname,double *qaold,ITG *mi,
             ITG *isolver,ITG *ncmat_,ITG *nstate_,ITG *iumat,
             double *cs,ITG *mcs,ITG *nkon,double **ener,ITG *mpcinfo,
             char *output,
             double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,
             double *physcon,ITG *nflow,double *ctrl,
             char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,ITG *nener,ITG *ikforc,ITG *ilforc,double *trab,
             ITG *inotr,ITG *ntrans,double **fmpcp,char *cbody,
             ITG *ibody,double *xbody,ITG *nbody,double *xbodyold,
             ITG *ielprop,double *prop,ITG *ntie,char *tieset,
	     ITG *itpamp,ITG *iviewfile,char *jobnamec,double *tietol,
	     ITG *nslavs,double *thicke,ITG *ics,
	     ITG *nintpoint,ITG *mortar,ITG *ifacecount,char *typeboun,
	     ITG **islavsurfp,double **pslavsurfp,double **clearinip,
	     ITG *nmat,double *xmodal,ITG *iaxial,ITG *inext);

void nonlingeo_precice(double **co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
	     ITG *ne,
	     ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	     ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,
             ITG *nmpc,
	     ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	     ITG **nelemloadp,char **sideloadp,double *xload,
	     ITG *nload,ITG *nactdof,
	     ITG **icolp,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
	     ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,ITG *ikboun,
	     ITG *ilboun,
	     double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	     double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
	     ITG **ielorienp,ITG *norien,double *orab,ITG *ntmat_,
	     double *t0,double *t1,double *t1old,
	     ITG *ithermal,double *prestr,ITG *iprestr,
	     double **vold,ITG *iperturb,double *sti,ITG *nzs, 
	     ITG *kode,char *filab,ITG *idrct,
	     ITG *jmax,ITG *jout,double *timepar,
	     double *eme,double *xbounold,
	     double *xforcold,double *xloadold,
             double *veold,double *accold,
             char *amname,double *amta,ITG *namta,ITG *nam,
             ITG *iamforc,ITG **iamloadp,
             ITG *iamt1,double *alpha,ITG *iexpl,
	     ITG *iamboun,double *plicon,ITG *nplicon,double *plkcon,
	     ITG *nplkcon,
             double **xstatep,ITG *npmat_,ITG *istep,double *ttime,
	     char *matname,double *qaold,ITG *mi,
             ITG *isolver,ITG *ncmat_,ITG *nstate_,ITG *iumat,
             double *cs,ITG *mcs,ITG *nkon,double **ener,ITG *mpcinfo,
             char *output,
             double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,
             double *physcon,ITG *nflow,double *ctrl,
             char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,ITG *nener,ITG *ikforc,ITG *ilforc,double *trab,
             ITG *inotr,ITG *ntrans,double **fmpcp,char *cbody,
             ITG *ibody,double *xbody,ITG *nbody,double *xbodyold,
             ITG *ielprop,double *prop,ITG *ntie,char *tieset,
	     ITG *itpamp,ITG *iviewfile,char *jobnamec,double *tietol,
	     ITG *nslavs,double *thicke,ITG *ics,
	     ITG *nintpoint,ITG *mortar,ITG *ifacecount,char *typeboun,
	     ITG **islavsurfp,double **pslavsurfp,double **clearinip,
	     ITG *nmat,double *xmodal,ITG *iaxial,ITG *inext,
		 char *preciceParticipantName,char *configFilename
					   );

void FORTRAN(nonlinmpc,(double *co,double *vold,ITG *ipompc,ITG *nodempc,
		   double *coefmpc,char *labmpc,ITG *nmpc,ITG *ikboun,
		   ITG *ilboun,ITG *nboun,double *xbounact,double *aux,
		   ITG *iaux,ITG *maxlenmpc,ITG *ikmpc,ITG *ilmpc,
                   ITG *icascade,ITG *kon,ITG *ipkon,char *lakon,
		   ITG *ne,double *reltime,ITG *newstep,double *xboun,
		   double *fmpc,ITG *newinc,ITG *idiscon,ITG *ncont,
		   double *trab,ITG *ntrans,ITG *ithermal,ITG *mi));

void FORTRAN(norm,(double *vel,double *velnorm,ITG *nef));

void FORTRAN(normalsoninterface,(ITG *istartset,ITG *iendset,
	     ITG *ialset,ITG *imast,ITG *ipkon,ITG *kon,char *lakon,
             ITG *imastnode,ITG *nmastnode,double *xmastnor,double *co));

void objectivesmmain_se(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	     ITG *ne,double *v,double *stn,ITG *inum,
	     double *stx,
	     double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	     double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	     ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	     double *t0,double *t1,ITG *ithermal,double *prestr,
             ITG *iprestr,char *filab,double *eme,double *emn,
             double *een,ITG *iperturb,double *f,double *fn,ITG *nactdof,
             ITG *iout,double *qa,
	     double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
	     double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,
	     double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod,
             double *cam,ITG *neq,double *veold,double *accold,
	     double *bet,double *gam,double *dtime,double *time,
             double *ttime,double *plicon,
             ITG *nplicon,double *plkcon,ITG *nplkcon,
             double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
	     double *epn,char *matname,ITG *mi,ITG *ielas,
	     ITG *icmd,ITG *ncmat_,ITG *nstate_,double *stiini,
	     double *vini,ITG *ikboun,ITG *ilboun,double *ener,
	     double *enern,double *emeini,double *xstaten,double *eei,
             double *enerini,double *cocon,ITG *ncocon,char *set,
             ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,double *qfx,double *qfn,double *trab,
             ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,
	     ITG *nload,ITG *ikmpc,ITG *ilmpc,ITG *istep,ITG *iinc,
	     double *springarea,double *reltime,ITG *ne0,double *xforc,
             ITG *nforc,double *thicke,
             double *shcon,ITG *nshcon,char *sideload,double *xload,
             double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
             double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
             ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
             ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
             double *energy,double *distmin,
	     ITG *ndesi,ITG *nodedesi,ITG *ndirdesi,ITG *nobject,
	     char *objectset,double *g0,double *dgdx,double *dgdv,double *sti,
	     double *dgdxtot,double *dfextminds,double *df);

void *objectivemt_shapeener(ITG *i);

void *objectivemt_mass(ITG *i);

void FORTRAN(op,(ITG *n,double *x,double *y,double *ad,double *au,ITG *jq,ITG *irow));

void FORTRAN(opas,(ITG *n,double *x,double *y,double *ad,double *au,ITG *jq,
		   ITG *irow,ITG *nzs));

void FORTRAN(op_corio,(ITG *n,double *x,double *y,double *ad,double *au,
		       ITG *jq,ITG *irow));

void FORTRAN(openfile,(char *jobname,char *output));

void FORTRAN(openfilefluid,(char *jobname));

void FORTRAN(postview,(ITG *ntr,char *sideload,ITG *nelemload,ITG *kontri,
             ITG *ntri,ITG *nloadtr,double *tenv,double *adview,double *auview,
             double *area,double *fenv,ITG *jqrad,ITG *irowrad,ITG *nzsrad));

void FORTRAN(precfd,(ITG *ne,ITG *ipkon,ITG *kon,char *lakon,ITG *ipnei,
                     ITG *neifa,ITG *neiel,ITG *ipoface,ITG *nodface,
                     ITG *ielfa,ITG *nkonnei,ITG *nface,ITG *ifaext,
                     ITG *nfaext,ITG *isolidsurf,ITG *nsolidsurf,char *set,
                     ITG *nset,ITG *istartset,ITG *iendset,ITG *ialset,
                     double *vel,double *vold,ITG *mi,ITG *neij,ITG *nef,
                     ITG *nactdoh,ITG *ipkonf,char *lakonf,ITG *ielmatf,
                     ITG *ielmat,ITG *ielorienf,ITG *ielorien,ITG *norien));

void precontact(ITG *ncont,ITG *ntie,char *tieset,ITG *nset,char *set,
        ITG *istartset,ITG *iendset,ITG *ialset,ITG *itietri,
        char *lakon,ITG *ipkon,ITG *kon,ITG *koncont,ITG *ne,
        double *cg,double *straight,double *co,double *vold,
        ITG *istep,ITG *iinc,ITG *iit,ITG *itiefac,
        ITG *islavsurf,ITG *islavnode,ITG *imastnode,
        ITG *nslavnode,ITG *nmastnode,ITG *imastop,ITG *mi,
	ITG *ipe,ITG *ime,double *tietol,ITG *iflagact,
	ITG *nintpoint,double **pslavsurfp,double *xmastnor,double *cs,
	ITG *mcs,ITG *ics,double *clearini,ITG *nslavs);

void prediction(double *uam,ITG *nmethod,double *bet,double *gam,double *dtime,
               ITG *ithermal,ITG *nk,double *veold,double *accold,double *v,
	       ITG *iinc,ITG *idiscon,double *vold,ITG *nactdof,ITG *mi);

void prediction_em(double *uam,ITG *nmethod,double *bet,double *gam,double *dtime,
               ITG *ithermal,ITG *nk,double *veold,double *v,
	       ITG *iinc,ITG *idiscon,double *vold,ITG *nactdof,ITG *mi);

void preiter(double *ad,double **aup,double *b,ITG **icolp,ITG **irowp,
	     ITG *neq,ITG *nzs,ITG *isolver,ITG *iperturb);

void FORTRAN(printout,(char *set,ITG *nset,ITG *istartset,ITG *iendset,
             ITG *ialset,ITG *nprint,char *prlab,char *prset,
             double *v,double *t1,double *fn,ITG *ipkon,
             char *lakon,double *stx,double *eei,double *xstate,
             double *ener,ITG *mi,ITG *nstate_,ITG *ithermal,
             double *co,ITG *kon,double *qfx,double *ttime,double *trab,
             ITG *inotr,ITG *ntrans,double *orab,ITG *ielorien,
	     ITG *norien,ITG *nk,ITG *ne,ITG *inum,char *filab,double *vold,
             ITG *ikin,ITG *ielmat,double *thicke,double *eme,ITG *islavsurf,
             ITG *mortar));

void FORTRAN(printoutfluid,(char *set,ITG *nset,ITG *istartset,ITG *iendset,
             ITG *ialset,ITG *nprint,char *prlab,char *prset,
             double *v,double *t1,double *fn,ITG *ipkon,
             char *lakon,double *stx,double *eei,double *xstate,
             double *ener,ITG *mi,ITG *nstate_,ITG *ithermal,
             double *co,ITG *kon,double *qfx,double *ttime,double *trab,
             ITG *inotr,ITG *ntrans,double *orab,ITG *ielorien,
	     ITG *norien,ITG *nk,ITG *ne,ITG *inum,char *filab,double *vold,
             ITG *ielmatf,double *thicke,double *eme,double *vcontu,
	     double *physcon,ITG *nactdoh,ITG *ielpropf,double *prop));

void FORTRAN(printoutface,(double *co,double *rhcon,ITG *nrhcon,ITG *ntmat_,
            double *vold,double *shcon,ITG *nshcon,double *cocon,
            ITG *ncocon,ITG *icompressible,ITG *istartset,ITG *iendset,
            ITG *ipkon,char *lakon,ITG *kon,ITG *ialset,char *prset,
	    double *timef,ITG *nset,char *set,ITG *nprint,char *prlab,
	    ITG *ielmat,ITG *mi));

int pthread_create (pthread_t *thread_id,const pthread_attr_t *attributes,
                    void *(*thread_function)(void *),void *arguments);

int pthread_join (pthread_t thread,void **status_ptr);

void radcyc(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
	    double *cs,ITG *mcs,ITG *nkon,ITG *ialset,ITG *istartset,
            ITG *iendset,ITG **kontrip,ITG *ntri,
            double **cop,double **voldp,ITG *ntrit,ITG *inocs,ITG *mi);

void radflowload(ITG *itg,ITG *ieg,ITG *ntg,ITG *ntr,double *adrad,
       double *aurad,
       double *bcr,ITG *ipivr,double *ac,double *bc,ITG *nload,
       char *sideload,ITG *nelemload,double *xloadact,char *lakon,ITG *ipiv,
       ITG *ntmat_,double *vold,double *shcon,ITG *nshcon,ITG *ipkon,
       ITG *kon,double *co,ITG *kontri,ITG *ntri,
       ITG *nloadtr,double *tarea,double *tenv,double *physcon,double *erad,
       double **adviewp,double **auviewp,
       ITG *nflow,ITG *ikboun,double *xboun,ITG *nboun,ITG *ithermal,
       ITG *iinc,ITG *iit,double *cs,ITG *mcs,ITG *inocs,ITG *ntrit,
       ITG *nk,double *fenv,ITG *istep,double *dtime,double *ttime,
       double *time,ITG *ilboun,ITG *ikforc,ITG *ilforc,double *xforc,
       ITG *nforc,double *cam,ITG *ielmat,ITG *nteq,double *prop,
       ITG *ielprop,ITG *nactdog,ITG *nacteq,ITG *nodeboun,ITG *ndirboun,
       ITG *network,double *rhcon,ITG *nrhcon,
       ITG *ipobody,ITG *ibody,double *xbody,ITG *nbody,ITG *iviewfile,
       char *jobnamef,double *ctrl,double *xloadold,double *reltime,
       ITG *nmethod,char *set,ITG *mi,ITG * istartset,ITG* iendset,
       ITG *ialset,ITG *nset,ITG *ineighe,ITG *nmpc,ITG *nodempc,
       ITG *ipompc,double *coefmpc,char *labmpc,ITG *iemchange,ITG *nam,
       ITG *iamload,ITG *jqrad,ITG *irowrad,ITG *nzsrad,ITG *icolrad,
       ITG *ne,ITG *iaxial,double *qa);

void FORTRAN (radmatrix,(ITG *ntr,double *adrad,double *aurad,double *bcr,
       char *sideload,ITG *nelemload,double *xloadact,char *lakon,
       double *vold,ITG *ipkon,ITG *kon,double *co,ITG *nloadtr,
       double *tarea,double *tenv,double *physcon,double *erad,
       double *adview,double *auview,ITG *ithermal,ITG *iinc,
       ITG *iit,double *fenv,ITG *istep,
       double *dtime,double *ttime,double *time,ITG *iviewfile,
       double *xloadold,double *reltime,ITG *nmethod,
       ITG *mi,ITG *iemchange,ITG *nam,ITG *iamload,ITG *jqrad,
       ITG *irowrad,ITG *nzsrad));

void FORTRAN (radresult,(ITG *ntr,double *xloadact,double *bcr,
       ITG *nloadtr,double *tarea,double * tenv,double *physcon,double *erad,
       double *auview,double *fenv,ITG *irowrad,ITG *jqrad,
       ITG *nzsrad,double *q));

void FORTRAN(readforce,(double *zc,ITG *neq,ITG *nev,ITG *nactdof,
	     ITG *ikmpc,ITG *nmpc,ITG *ipompc,ITG *nodempc,ITG *mi,
	     double *coefmpc,char *jobnamec,double *aa,
             ITG *igeneralizedforce));

void readinput(char *jobnamec,char **inpcp,ITG *nline,ITG *nset,ITG *ipoinp,
        ITG **inpp,ITG **ipoinpcp,ITG *ithermal); 

void FORTRAN(readview,(ITG *ntr,double *adview,double *auview,double *fenv,
             ITG *nzsrad,ITG *ithermal,char *jobnamef));

void FORTRAN(rearrange,(double *au,ITG *irow,ITG *icol,ITG *ndim,ITG *neq));

void FORTRAN(rectcyl,(double *co,double *v,double *fn,double *stn,
		      double *qfn,double *een,double *cs,ITG *nk,
                      ITG *icntrl,double *t,char *filab,ITG *imag,
                      ITG *mi,double *emn));

void FORTRAN(rectcylexp,(double *co,double *v,double *fn,double *stn,
		      double *qfn,double *een,double *cs,ITG *nkt,
		      ITG *icntrl,double *t,char *filab,ITG *imag,ITG *mi,
		      ITG *iznode,ITG *nznode,ITG *nsectors,ITG *nk,
                      double *emn));

void FORTRAN(rectcyltrfm,(ITG *node,double *co,double *cs,ITG *cntrl,
             double *fin,double *fout));

void FORTRAN(rectcylvi,(double *co,double *v,double *fn,double *stn,
		      double *qfn,double *een,double *cs,ITG *nk,
		      ITG *icntrl,double *t,char *filab,ITG *imag,ITG *mi,
                      double *emn));

void remastruct(ITG *ipompc,double **coefmpcp,ITG **nodempcp,ITG *nmpc,
              ITG *mpcfree,ITG *nodeboun,ITG *ndirboun,ITG *nboun,
              ITG *ikmpc,ITG *ilmpc,ITG *ikboun,ITG *ilboun,
              char *labmpc,ITG *nk,
              ITG *memmpc_,ITG *icascade,ITG *maxlenmpc,
              ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
              ITG *nactdof,ITG *icol,ITG *jq,ITG **irowp,ITG *isolver,
              ITG *neq,ITG *nzs,ITG *nmethod,double **fp,
              double **fextp,double **bp,double **aux2p,double **finip,
              double **fextinip,double **adbp,double **aubp,ITG *ithermal,
	      ITG *iperturb,ITG *mass,ITG *mi,ITG *iexpl,ITG *mortar,
	      char *typeboun,double **cvp,double **cvinip,ITG *iit);

void remastructar(ITG *ipompc,double **coefmpcp,ITG **nodempcp,ITG *nmpc,
              ITG *mpcfree,ITG *nodeboun,ITG *ndirboun,ITG *nboun,
              ITG *ikmpc,ITG *ilmpc,ITG *ikboun,ITG *ilboun,
              char *labmpc,ITG *nk,
              ITG *memmpc_,ITG *icascade,ITG *maxlenmpc,
              ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
              ITG *nactdof,ITG *icol,ITG *jq,ITG **irowp,ITG *isolver,
              ITG *neq,ITG *nzs,ITG *nmethod,ITG *ithermal,
	      ITG *iperturb,ITG *mass,ITG *mi,ITG *ics,double *cs,
	      ITG *mcs,ITG *mortar,char *typeboun,ITG *iit);

void remastructem(ITG *ipompc,double **coefmpcp,ITG **nodempcp,ITG *nmpc,
              ITG *mpcfree,ITG *nodeboun,ITG *ndirboun,ITG *nboun,
              ITG *ikmpc,ITG *ilmpc,ITG *ikboun,ITG *ilboun,
              char *labmpc,ITG *nk,
              ITG *memmpc_,ITG *icascade,ITG *maxlenmpc,
              ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
              ITG *nactdof,ITG *icol,ITG *jq,ITG **irowp,ITG *isolver,
              ITG *neq,ITG *nzs,ITG *nmethod,double **fp,
              double **fextp,double **bp,double **aux2p,double **finip,
              double **fextinip,double **adbp,double **aubp,ITG *ithermal,
	      ITG *iperturb,ITG *mass,ITG *mi,ITG *ielmat,double *elcon,
	      ITG *ncmat_,ITG *ntmat_,ITG *inomat);

void FORTRAN(restartshort,(ITG *nset,ITG *nload,ITG *nbody,ITG *nforc,
    ITG *nboun,
    ITG *nk,ITG *ne,ITG *nmpc,ITG *nalset,ITG *nmat,ITG *ntmat,ITG *npmat,
    ITG *norien,ITG *nam,ITG *nprint,ITG *mint,ITG *ntrans,ITG *ncs,
    ITG *namtot,ITG *ncmat,ITG *memmpc,ITG *ne1d,ITG *ne2d,ITG *nflow,
    char *set,ITG *meminset,ITG *rmeminset,char *jobnamec,ITG *irestartstep,
    ITG *icntrl,ITG *ithermal,ITG *nener,ITG *nstate_,ITG *ntie,ITG *nslavs,
    ITG *nkon,ITG *mcs,ITG *nprop,ITG *mortar,ITG *ifacecount,ITG *nintpoint,
    ITG *infree));

void FORTRAN(restartwrite,(ITG *istep,ITG *nset,ITG*nload,ITG *nforc,
  ITG * nboun,ITG *nk,ITG *ne,ITG *nmpc,ITG *nalset,ITG *nmat,ITG *ntmat_,
  ITG *npmat_,ITG *norien,ITG *nam,ITG *nprint,ITG *mi,
  ITG *ntrans,ITG *ncs_,ITG *namtot_,ITG *ncmat_,ITG *mpcend,
  ITG *maxlenmpc,ITG *ne1d,
  ITG *ne2d,ITG *nflow,ITG *nlabel,ITG *iplas,ITG *nkon,ITG *ithermal,
  ITG *nmethod,ITG *iperturb,ITG *nstate_,ITG *nener,char *set,
  ITG *istartset,ITG *iendset,ITG *ialset,double *co,ITG *kon,ITG *ipkon,
  char *lakon,ITG *nodeboun,ITG *ndirboun,ITG *iamboun,double *xboun,
  ITG *ikboun,ITG *ilboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
  char *labmpc,ITG *ikmpc,ITG *ilmpc,ITG *nodeforc,ITG *ndirforc,
  ITG *iamforc,double *xforc,ITG *ikforc,ITG *ilforc,ITG *nelemload,
  ITG *iamload,char *sideload,double *xload, 
  double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,double *alcon,
  ITG *nalcon,double *alzero,double *plicon,ITG *nplicon,double *plkcon,
  ITG *nplkcon,char *orname,double *orab,ITG *ielorien,double *trab,
  ITG *inotr,char *amname,double *amta,ITG *namta,double *t0,double *t1,
  ITG *iamt1,double *veold,ITG *ielmat,char *matname,
  char *prlab,char *prset,char *filab,double *vold,
  ITG *nodebounold,ITG *ndirbounold,double *xbounold,double *xforcold,
  double *xloadold,double *t1old,double *eme,ITG *iponor,
  double *xnor,ITG *knor,double *thicke,double *offset,
  ITG *iponoel,ITG *inoel,ITG *rig,
  double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,
  ITG *ics,double *sti,double *ener,double *xstate,
  char *jobnamec,ITG *infree,double *prestr,ITG *iprestr,
  char *cbody,ITG *ibody,double *xbody,ITG *nbody,double *xbodyold,
  double *ttime,double *qaold,double *cs,
  ITG *mcs,char *output,double *physcon,double *ctrl,char *typeboun,
  double *fmpc,char *tieset,ITG *ntie,double *tietol,ITG *nslavs,
  double *t0g,double *t1g,ITG *nprop,ITG *ielprop,double *prop,ITG *mortar,
  ITG *nintpoint,ITG *ifacecount,ITG *islavsurf,double *pslavsurf,
  double *clearini));

void FORTRAN(resultnet,(ITG *itg,ITG *ieg,ITG *ntg,
                        double *bc,ITG *nload,char *sideload,
                        ITG *nelemload,double *xloadact,char *lakon,
                        ITG *ntmat_,double *v,double *shcon,ITG *nshcon,
                        ITG *ipkon,ITG *kon,double *co,ITG *nflow,
			ITG *iinc,ITG *istep,
                        double *dtime,double *ttime,double *time,
			ITG *ikforc,ITG *ilforc,
                        double *xforcact,ITG *nforc,
                        ITG *ielmat,ITG *nteq,double *prop,
                        ITG *ielprop,ITG *nactdog,ITG *nacteq,ITG *iin,
                        double *physcon,double *camt,double *camf,
                        double *camp,double *rhcon,ITG *nrhcon,
			ITG *ipobody,ITG *ibody,double *xbody,ITG *nbody,
                        double *dtheta,double *vold,double *xloadold,
                        double *reltime,ITG *nmethod,char *set,ITG *mi,
                        ITG *ineighe,double *cama,double *vamt,
                        double *vamf,double *vamp,double *vama,
                        ITG *nmpc,ITG *nodempc,ITG *ipompc,double *coefmpc,
                        char *labmpc,ITG *iaxial));

void results(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	     ITG *ne,double *v,double *stn,ITG *inum,
	     double *stx,
	     double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	     double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	     ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	     double *t0,double *t1,ITG *ithermal,double *prestr,
             ITG *iprestr,char *filab,double *eme,double *emn,
             double *een,ITG *iperturb,double *f,double *fn,ITG *nactdof,
             ITG *iout,double *qa,
	     double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
	     double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,
	     double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod,
             double *vmax,ITG *neq,double *veold,double *accold,
	     double *beta,double *gamma,double *dtime,double *time,
             double *ttime,double *plicon,
             ITG *nplicon,double *plkcon,ITG *nplkcon,
             double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
	     double *epl,char *matname,ITG *mi,ITG *ielas,
	     ITG *icmd,ITG *ncmat_,ITG *nstate_,double *stiini,
	     double *vini,ITG *ikboun,ITG *ilboun,double *ener,
	     double *enern,double *emeini,double *xstaten,double *eei,
             double *enerini,double *cocon,ITG *ncocon,char *set,
             ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,double *qfx,double *qfn,double *trab,
             ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,
	     ITG *nload,ITG *ikmpc,ITG *ilmpc,ITG *istep,ITG *iinc,
	     double *springarea,double *reltime,ITG *ne0,double *xforc,
             ITG *nforc,double *thicke,
             double *shcon,ITG *nshcon,char *sideload,double *xload,
             double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
             double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
             ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
             ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
             double *energy,ITG *kscale);

void FORTRAN(resultsem,(double *co,ITG *kon,ITG *ipkon,char *lakon,
             double *v,double *elcon,ITG *nelcon,ITG *ielmat,ITG *ntmat_,
             double *vold,double *dtime,char *matname,ITG *mi,ITG *ncmat_,
             ITG *nea,ITG *neb,double *sti,double *alcon,
	     ITG *nalcon,double *h0,ITG *istartset,ITG *iendset,ITG *ialset,
	     ITG *iactive,double *fn));

void *resultsemmt(ITG *i);

void  FORTRAN(resultsforc,(ITG *nk,double *f,double *fn,ITG *nactdof,
       ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,ITG *nmpc,
       ITG *mi,double *fmpc,ITG *calcul_fn,
       ITG *calcul_f));

void  FORTRAN(resultsforc_em,(ITG *nk,double *f,double *fn,ITG *nactdof,
       ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,ITG *nmpc,
       ITG *mi,double *fmpc,ITG *calcul_fn,ITG *calcul_f,ITG *inomat));
       
void  FORTRAN(resultsforc_se,(ITG *nk,double *f,double *dfn,ITG *nactdof,
       ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,ITG *nmpc,
       ITG *mi,double *fmpc,ITG *calcul_fn,ITG *calcul_f,ITG *ndesi,
       double *df));

void FORTRAN(resultsini,(ITG *nk,double *v,ITG *ithermal,char *filab,
       ITG *iperturb,double *f,double *fn,ITG *nactdof,ITG *iout,
       double *qa,double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
       double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
       char *labmpc,ITG *nmpc,ITG *nmethod,double *cam,ITG *neq,
       double *veold,double *accold,double *bet,double *gam,double *dtime,
       ITG *mi,double *vini,ITG *nprint,char *prlab,ITG *intpointvar,
       ITG *calcul_fn,ITG *calcul_f,ITG *calcul_qa,ITG *calcul_cauchy,
       ITG *iener,ITG *ikin,ITG *intpointvart,double *xforc,ITG *nforc));

void FORTRAN(resultsini_em,(ITG *nk,double *v,ITG *ithermal,char *filab,
       ITG *iperturb,double *f,double *fn,ITG *nactdof,ITG *iout,
       double *qa,double *b,ITG *nodeboun,ITG *ndirboun,
       double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
       char *labmpc,ITG *nmpc,ITG *nmethod,double *cam,ITG *neq,
       double *veold,double *dtime,
       ITG *mi,double *vini,ITG *nprint,char *prlab,ITG *intpointvar,
       ITG *calcul_fn,ITG *calcul_f,ITG *calcul_qa,ITG *calcul_cauchy,
       ITG *iener,ITG *ikin,ITG *intpointvart,double *xforc,ITG *nforc));

void FORTRAN(resultsmech,(double *co,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
          double *v,double *stx,double *elcon,ITG *nelcon,double *rhcon,
          ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
          ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,double *t0,
          double *t1,ITG *ithermal,double *prestr,ITG *iprestr,double *eme,
          ITG *iperturb,double *fn,ITG *iout,double *qa,double *vold,
          ITG *nmethod,double *veold,double *dtime,double *time,
          double *ttime,double *plicon,ITG *nplicon,double *plkcon,
          ITG *nplkcon,double *xstateini,double *xstiff,double *xstate,
          ITG *npmat_,char *matname,ITG *mi,ITG *ielas,ITG *icmd,ITG *ncmat_,
          ITG *nstate_,double *stiini,double *vini,double *ener,double *eei,
          double *enerini,ITG *istep,ITG *iinc,double *springarea,
          double *reltime,ITG *calcul_fn,ITG *calcul_qa,ITG *calcul_cauchy,
	  ITG *iener,ITG *ikin,ITG *nal,ITG *ne0,double *thicke,
	  double *emeini,double *pslavsurf,double *pmastsurf,ITG *mortar,
	  double *clearini,ITG *nea,ITG *neb,ITG *ielprop,double *prop,
          ITG *kscale));

void *resultsmechmt(ITG *i);

void *resultsmechmt_se(ITG *i);

void FORTRAN(resultsmech_se,(double *co,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
          double *v,double *stx,double *elcon,ITG *nelcon,double *rhcon,
          ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
          ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,double *t0,
          double *t1,ITG *ithermal,double *prestr,ITG *iprestr,double *eme,
          ITG *iperturb,double *fn,ITG *iout,double *qa,double *vold,
          ITG *nmethod,double *veold,double *dtime,double *time,
          double *ttime,double *plicon,ITG *nplicon,double *plkcon,
          ITG *nplkcon,double *xstateini,double *xstiff,double *xstate,
          ITG *npmat_,char *matname,ITG *mi,ITG *ielas,ITG *icmd,ITG *ncmat_,
          ITG *nstate_,double *stiini,double *vini,double *ener,double *eei,
          double *enerini,ITG *istep,ITG *iinc,double *springarea,
          double *reltime,ITG *calcul_fn,ITG *calcul_qa,ITG *calcul_cauchy,
	  ITG *iener,ITG *ikin,ITG *nal,ITG *ne0,double *thicke,
	  double *emeini,double *pslavsurf,double *pmastsurf,ITG *mortar,
          double *clearini,ITG *nea,ITG *neb,ITG *ielprop,double *prop,
	  double *dfn,double *distmin,ITG *ndesi,ITG *nodedesi,
	  ITG *ndirdesi,double *fn0,double *sti));

void  FORTRAN(resultsprint,(double *co,ITG *nk,ITG *kon,ITG *ipkon,
       char *lakon,ITG *ne,double *v,double *stn,ITG *inum,double *stx,
       ITG *ielorien,ITG *norien,double *orab,double *t1,ITG *ithermal,
       char *filab,double *een,ITG *iperturb,double *fn,ITG *nactdof,
       ITG *iout,double *vold,ITG *nodeboun,ITG *ndirboun,ITG *nboun,
       ITG *nmethod,double *ttime,double *xstate,double *epn,ITG *mi,
       ITG *nstate_,double *ener,double *enern,double *xstaten,double *eei,
       char *set,ITG *nset,ITG *istartset,ITG *iendset,ITG *ialset,ITG *nprint,
       char *prlab,char *prset,double *qfx,double *qfn,double *trab,ITG *inotr,
       ITG *ntrans,ITG *nelemload,ITG *nload,ITG *ikin,ITG *ielmat,
       double *thicke,double *eme,double *emn,double *rhcon,ITG *nrhcon,
       double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,ITG *ntmat_,
       char *sideload,ITG *icfd,ITG *inomat,double *pslavsurf,
       ITG *islavact,double *cdn,ITG *mortar,ITG *islavnode,ITG *nslavnode,
       ITG *ntie,ITG *islavsurf,double *time,ITG *ielprop,double *prop,
       double *veold,ITG *ne0,ITG *nmpc,ITG *ipompc,ITG *nodempc,
       char *labmpc,double *energyini,double *energy));
	     
void results_se(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	     ITG *ne,double *v,double *stn,ITG *inum,
	     double *stx,
	     double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	     double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	     ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	     double *t0,double *t1,ITG *ithermal,double *prestr,
             ITG *iprestr,char *filab,double *eme,double *emn,
             double *een,ITG *iperturb,double *f,double *fn,ITG *nactdof,
             ITG *iout,double *qa,
	     double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
	     double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,
	     double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod,
             double *vmax,ITG *neq,double *veold,double *accold,
	     double *beta,double *gamma,double *dtime,double *time,
             double *ttime,double *plicon,
             ITG *nplicon,double *plkcon,ITG *nplkcon,
             double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
	     double *epl,char *matname,ITG *mi,ITG *ielas,
	     ITG *icmd,ITG *ncmat_,ITG *nstate_,double *stiini,
	     double *vini,ITG *ikboun,ITG *ilboun,double *ener,
	     double *enern,double *emeini,double *xstaten,double *eei,
             double *enerini,double *cocon,ITG *ncocon,char *set,
             ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,double *qfx,double *qfn,double *trab,
             ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,
	     ITG *nload,ITG *ikmpc,ITG *ilmpc,ITG *istep,ITG *iinc,
	     double *springarea,double *reltime,ITG *ne0,double *xforc,
             ITG *nforc,double *thicke,
             double *shcon,ITG *nshcon,char *sideload,double *xload,
             double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
             double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
             ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
             ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
             double *energy,double *df,double *distmin,
	     ITG *ndesi,ITG *nodedesi,ITG *ndirdesi,double *sti);

void FORTRAN(resultstherm,(double *co,ITG *kon,ITG *ipkon,
       char *lakon,double *v,
       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,ITG *ielmat,
       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,double *t0,
       ITG *iperturb,double *fn,double *shcon,ITG *nshcon,ITG *iout,
       double *qa,double *vold,ITG *ipompc,ITG *nodempc,
       double *coefmpc,ITG *nmpc,double *dtime,
       double *time,double *ttime,double *plkcon,ITG *nplkcon,double *xstateini,
       double *xstiff,double *xstate,ITG *npmat_,char *matname,
       ITG *mi,ITG *ncmat_,ITG *nstate_,double *cocon,ITG *ncocon,
       double *qfx,ITG *ikmpc,ITG *ilmpc,ITG *istep,
       ITG *iinc,double *springarea,ITG *calcul_fn,ITG *calcul_qa,ITG *nal,
       ITG *nea,ITG *neb,ITG *ithermal,ITG *nelemload,ITG *nload,
       ITG *nmethod,double *reltime,char *sideload,double *xload,
       double *xloadold,double *pslavsurf,double *pmastsurf,ITG *mortar,
       double *clearini,double *plicon,ITG *nplicon,ITG *ielprop,double *prop));

void *resultsthermemmt(ITG *i);

void *resultsthermmt(ITG *i);

void *resultsthermmt_se(ITG *i);

void resultsinduction(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	     ITG *ne,double *v,double *stn,ITG *inum,
	     double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	     double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	     ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	     double *t0,double *t1,ITG *ithermal,double *prestr,
             ITG *iprestr,char *filab,double *eme,double *emn,
             double *een,ITG *iperturb,double *f,double *fn,ITG *nactdof,
             ITG *iout,double *qa,
	     double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
	     double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,
	     double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod,
             double *vmax,ITG *neq,double *veold,double *accold,
	     double *beta,double *gamma,double *dtime,double *time,
             double *ttime,double *plicon,
             ITG *nplicon,double *plkcon,ITG *nplkcon,
             double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
	     double *epl,char *matname,ITG *mi,ITG *ielas,
	     ITG *icmd,ITG *ncmat_,ITG *nstate_,double *sti,
	     double *vini,ITG *ikboun,ITG *ilboun,double *ener,
	     double *enern,double *emeini,double *xstaten,double *eei,
             double *enerini,double *cocon,ITG *ncocon,char *set,
             ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,double *qfx,double *qfn,double *trab,
             ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,
	     ITG *nload,ITG *ikmpc,ITG *ilmpc,ITG *istep,ITG *iinc,
	     double *springarea,double *reltime,ITG *ne0,double *xforc,
             ITG *nforc,double *thicke,
             double *shcon,ITG *nshcon,char *sideload,double *xload,
	     double *xloadold,ITG *icfd,ITG *inomat,double *h0,
	     ITG *islavnode,ITG *nslavnode,ITG *ntie,ITG *ielprop,
	     double *prop,ITG *iactive,double *energy,double *energyini);

void FORTRAN(rhs,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,ITG *ipompc,ITG *nodempc,double *coefmpc,
	       ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,double *bb,ITG *nactdof,ITG *neq,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,
               ITG *iprestr,double *vold,ITG *iperturb,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               ITG *npmat_,double *ttime,double *time,ITG *istep,
               ITG *iinc,double *dtime,double *physcon,ITG *ibody,
	       double *xbodyold,double *reltime,double *veold,
	       char *matname,ITG *mi,ITG *ikactmech,ITG *nactmech,
               ITG *ielprop,double *prop));

void FORTRAN(rhsp,(ITG *ne,char *lakon,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
	     double *xle,double *b,double *xxn,
	     ITG *neq,ITG *nzs,double *hfa,double *gradpel,double *bp,
	     double *xxi,ITG *neij,double *xlen,ITG *nefa,ITG *nefb));

void rhspmain(ITG *ne,char *lakon,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
	     double *xle,double *b,double *xxn,
	     ITG *neq,ITG *nzs,double *hfa,double *gradpel,double *bp,
	     double *xxi,ITG *neij,double *xlen,ITG *iatleastonepressurebc);

void *rhspmt(ITG *i);
       
void sensitivity(double *co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
	     ITG *ne,
	     ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	     ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
             ITG *nmpc,
	     ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	     ITG *nelemload,char *sideload,double *xload,
	     ITG *nload,ITG *nactdof,
	     ITG **icolp,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
	     ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	     ITG *ilboun,
	     double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	     double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
	     ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	     double *t0,double *t1,double *t1old,
	     ITG *ithermal,double *prestr,ITG *iprestr,
	     double *vold,ITG *iperturb,double *sti,ITG *nzs,
	     ITG *kode,char *filab,double *eme,
             ITG *iexpl,double *plicon,ITG *nplicon,double *plkcon,
             ITG *nplkcon,
             double **xstatep,ITG *npmat_,char *matname,ITG *isolver,
	     ITG *mi,ITG *ncmat_,ITG *nstate_,double *cs,ITG *mcs,
             ITG *nkon,double **enerp,double *xbounold,
	     double *xforcold,double *xloadold,
             char *amname,double *amta,ITG *namta,
             ITG *nam,ITG *iamforc,ITG *iamload,
             ITG *iamt1,ITG *iamboun,double *ttime,char *output,
             char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,ITG *nener,double *trab,
             ITG *inotr,ITG *ntrans,double *fmpc,char *cbody,ITG *ibody,
	     double *xbody,ITG *nbody,double *xbodyold,double *timepar,
	     double *thicke,char *jobnamec,char *tieset,ITG *ntie,
	     ITG *istep,ITG *nmat,ITG *ielprop,double *prop,char *typeboun,
	     ITG *mortar,ITG *mpcinfo,double *tietol,ITG *ics,
	     ITG *icontact,ITG *nobject,char *objectset,ITG *istat);

void sensitivity_glob(double *dgdxtot,double *dgdxtotglob,ITG *nobject,
                  ITG *ndesi,ITG *nodedesi,ITG *ndirdesi,ITG *nk);
		  
void FORTRAN(shape3tri,(double *xi,double *et,double *xl,double *xsj,
                      double *xs,double *shp,ITG *iflag));

void FORTRAN(shape4q,(double *xi,double *et,double *xl,double *xsj,
                      double *xs,double *shp,ITG *iflag));

void FORTRAN(shape4tet,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,ITG *iflag));

void FORTRAN(shape6tri,(double *xi,double *et,double *xl,double *xsj,
                      double *xs,double *shp,ITG *iflag));

void FORTRAN(shape6w,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,ITG *iflag));

void FORTRAN(shape8h,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,ITG *iflag));

void FORTRAN(shape8q,(double *xi,double *et,double *xl,double *xsj,
                      double *xs,double *shp,ITG *iflag));

void FORTRAN(shape10tet,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,ITG *iflag));

void FORTRAN(shape15w,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,ITG *iflag));

void FORTRAN(shape20h,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,ITG *iflag));
       
void FORTRAN(slavintpoints,(ITG *ntie,ITG *itietri,ITG *ipkon,
        ITG *kon,char *lakon,double *straight,
        ITG *nintpoint,ITG *koncont,double *co,double *vold,double *xo,
        double *yo,double *zo,double *x,double *y,double *z,ITG *nx,
        ITG *ny,ITG *nz,ITG *islavsurf,
        ITG *islavnode,ITG *nslavnode,ITG *imastop,
	ITG *mi,ITG *ncont,ITG *ipe,ITG *ime,double *pslavsurf,
        ITG *i,ITG *l,ITG *ntri));

void FORTRAN(sortev,(ITG *nev,ITG *nmd,double *eigxx,ITG *cyclicsymmetry,
		     double *xx,double *eigxr,ITG *pev,
		     ITG *istartnmd,ITG *iendnmd,double *aa,double *bb));

void FORTRAN(spcmatch,(double *xboun,ITG *nodeboun,ITG *ndirboun,ITG *nboun,
	       double *xbounold,ITG *nodebounold,ITG *ndirbounold,
	       ITG *nbounold,ITG *ikboun,ITG *ilboun,double *vold,
	       double *reorder,ITG *nreorder,ITG *mi));

void FORTRAN(splitline,(char *text,char *textpart,ITG *n));

void spooles(double *ad,double *au,double *adb,double *aub,
             double *sigma,double *b,
	     ITG *icol,ITG *irow,ITG *neq,ITG *nzs,ITG *symmtryflag,
             ITG *inputformat,ITG *nzs3);

void FORTRAN(springforc_n2f,(double *xl,ITG *konl,double *vl,ITG *imat,
             double *elcon,ITG *nelcon,double *elas,double *fnl,ITG *ncmat_,
             ITG *ntmat_,ITG *nope,char *lakonl,double *t1l,ITG *kode,
             double *elconloc,double *plicon,ITG *nplicon,ITG *npmat_,
             double *senergy,ITG *iener,double *cstr,ITG *mi,
             double *springarea,ITG *nmethod,ITG *ne0,ITG *nstate_,
	     double *xstateini,double *xstate,double *reltime,
	     ITG *ielas,double *venergy));

void FORTRAN(springstiff_n2f,(double *xl,double *elas,ITG *konl,double *voldl,
             double *s,ITG *imat,double *elcon,ITG *nelcon,ITG *ncmat_,
             ITG *ntmat_,ITG *nope,char *lakonl,double *t1l,ITG *kode,
             double *elconloc,double *plicon,ITG *nplicon,ITG *npmat_,
             ITG *iperturb,double *springarea,ITG *nmethod,ITG *mi,ITG *ne0,
             ITG *nstate_,double *xstateini,double *xstate,double *reltime,
             ITG *nasym));

void steadystate(double **co,ITG *nk,ITG **kon,ITG **ipkon,char **lakon,ITG *ne,
	  ITG **nodeboun,ITG **ndirboun,double **xboun,ITG *nboun,
	  ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,ITG *nmpc,
	  ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	  ITG *nelemload,char *sideload,double *xload,
	  ITG *nload,
	  ITG **nactdof,ITG *neq,ITG *nzl,ITG *icol,ITG *irow,
	  ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,ITG **ikboun,
	  ITG **ilboun,
	  double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
          double *cocon,ITG *ncocon,
	  double *alcon,ITG *nalcon,double *alzero,ITG **ielmat,
	  ITG **ielorien,ITG *norien,double *orab,ITG *ntmat_,
	  double **t0,
	  double **t1,ITG *ithermal,double *prestr,ITG *iprestr,
	  double **voldp,ITG *iperturb,double *sti,ITG *nzs,
	  double *timepar,double *xmodal,
	  double **veoldp,char *amname,double *amta,
	  ITG *namta,ITG *nam,ITG *iamforc,ITG *iamload,
	  ITG **iamt1,ITG *jout,ITG *kode,char *filab,
	  double **emep,double *xforcold,double *xloadold,
          double **t1old,ITG **iamboun,
          double **xbounold,ITG *iexpl,double *plicon,ITG *nplicon,
          double *plkcon,ITG *nplkcon,
          double *xstate,ITG *npmat_,char *matname,ITG *mi,
          ITG *ncmat_,ITG *nstate_,double **enerp,char *jobnamec,
          double *ttime,char *set,ITG *nset,ITG *istartset,
          ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
          char *prset,ITG *nener,double *trab,
          ITG **inotr,ITG *ntrans,double **fmpcp,char *cbody,ITG *ibody,
          double *xbody,ITG *nbody,double *xbodyold,ITG *istep,
          ITG *isolver,ITG *jq,char *output,ITG *mcs,ITG *nkon,
	  ITG *ics,double *cs,ITG *mpcend,double *ctrl,
	  ITG *ikforc,ITG *ilforc,double *thicke,ITG *nmat,
          char *typeboun,ITG *ielprop,double *prop);

void FORTRAN(stop,());

void storecontactdof(ITG *nope,ITG *nactdof,ITG *mt,ITG *konl,
          ITG **ikactcontp,
          ITG *nactcont,ITG *nactcont_,double *bcont,double *fnl,
          ITG *ikmpc,ITG *nmpc,ITG *ilmpc,ITG *ipompc,ITG *nodempc,
	  double *coefmpc);

void FORTRAN(storeresidual,(ITG *nactdof,double *b,double *fn,char *filab,
             ITG *ithermal,ITG *nk,double *sti,double *stn,
             ITG *ipkon,ITG *inum,ITG *kon,char *lakon,
             ITG *ne,ITG *mi,double *orab,ITG *ielorien,
             double *co,ITG *itg,ITG *ntg,double *vold,
	     ITG *ielmat,double *thicke,ITG *ielprop,double *prop));
     
void FORTRAN(storecontactprop,(ITG *ne,ITG *ne0,char *lakon,ITG *kon,
	 ITG *ipkon,ITG *mi,ITG *ielmat,double *elcon,ITG *mortar,
	 double *adb,ITG *nactdof,double *springarea,ITG *ncmat_,
         ITG *ntmat_,double *stx,double *temax));

ITG strcmp1(const char *s1,const char *s2);

ITG strcmp2(const char *s1,const char *s2,ITG length);

ITG strcpy1(char *s1,const char *s2,ITG length);

void FORTRAN(subspace,(double *d,double *aa,double *bb,double *cc,
             double *alpham,double *betam,ITG *nev,
             double *xini,double *cd,double *cv,double *time,
             double *rwork,ITG *lrw,ITG *k,ITG *jout,double *rpar,
	     double *bj,ITG *iwork,ITG *liw,ITG *iddebdf,double *bjp));

void FORTRAN(tempload,(double *xforcold,double *xforc,double *xforcact,
               ITG *iamforc,ITG *nforc,double *xloadold,double *xload,
               double *xloadact,ITG *iamload,ITG *nload,ITG *ibody,
               double *xbody,ITG *nbody,double *xbodyold,double *xbodyact,
               double *t1old,double *t1,double *t1act,ITG *iamt1,
               ITG *nk,double *amta,ITG *namta,ITG *nam,double *ampli,
               double *time,double *reltime,double *ttime,double *dtime,
               ITG *ithermal,ITG *nmethod,
	       double *xbounold,double *xboun,double *xbounact,
	       ITG *iamboun,ITG *nboun,ITG *nodeboun,
               ITG *ndirboun,ITG *nodeforc,ITG *ndirforc,ITG *istep,
               ITG *iint,double *co,double *vold,ITG *itg,ITG *ntg,
               char *amname,ITG *ikboun,ITG *ilboun,ITG *nelemload,
	       char *sideload,ITG *mi,ITG *ntrans,double *trab,
               ITG *inotr,double *veold,ITG *integerglob,
               double *doubleglob,char *tieset,ITG *istartset,
               ITG *iendset,ITG *ialset,ITG *ntie,ITG *nmpc,ITG *ipompc,
               ITG *ikmpc,ITG *ilmpc,ITG *nodempc,double *coefmpc));

void FORTRAN(tempload_em,(double *xforcold,double *xforc,double *xforcact,
               ITG *iamforc,ITG *nforc,double *xloadold,double *xload,
               double *xloadact,ITG *iamload,ITG *nload,ITG *ibody,
               double *xbody,ITG *nbody,double *xbodyold,double *xbodyact,
               double *t1old,double *t1,double *t1act,ITG *iamt1,
               ITG *nk,double *amta,ITG *namta,ITG *nam,double *ampli,
               double *time,double *reltime,double *ttime,double *dtime,
               ITG *ithermal,ITG *nmethod,
	       double *xbounold,double *xboun,double *xbounact,
	       ITG *iamboun,ITG *nboun,ITG *nodeboun,
               ITG *ndirboun,ITG *nodeforc,ITG *ndirforc,ITG *istep,
               ITG *iint,double *co,double *vold,ITG *itg,ITG *ntg,
               char *amname,ITG *ikboun,ITG *ilboun,ITG *nelemload,
	       char *sideload,ITG *mi,ITG *ntrans,double *trab,
               ITG *inotr,double *veold,ITG *integerglob,
               double *doubleglob,char *tieset,ITG *istartset,
               ITG *iendset,ITG *ialset,ITG *ntie,ITG *nmpc,ITG *ipompc,
	       ITG *ikmpc,ITG *ilmpc,ITG *nodempc,double *coefmpc,
	       double *h0scale,ITG *inomat));

void FORTRAN(temploaddiff,(double *xforcold,double *xforc,double *xforcact,
               ITG *iamforc,ITG *nforc,double *xloadold,double *xload,
               double *xloadact,ITG *iamload,ITG *nload,ITG *ibody,
               double *xbody,ITG *nbody,double *xbodyold,double *xbodyact,
               double *t1old,double *t1,double *t1act,ITG *iamt1,
               ITG *nk,double *amta,ITG *namta,ITG *nam,double *ampli,
               double *time,double *reltime,double *ttime,double *dtime,
               ITG *ithermal,ITG *nmethod,
	       double *xbounold,double *xboun,double *xbounact,
	       ITG *iamboun,ITG *nboun,ITG *nodeboun,
               ITG *ndirboun,ITG *nodeforc,ITG *ndirforc,ITG *istep,
               ITG *iint,double *co,double *vold,ITG *itg,ITG *ntg,
               char *amname,ITG *ikboun,ITG *ilboun,ITG *nelemload,
	       char *sideload,ITG *mi,double *xforcdiff,double *xloaddiff,
	       double *xbodydiff,double *t1diff,double *xboundiff,
	       ITG *icorrect,ITG *iprescribedboundary,ITG *ntrans,
               double *trab,ITG *inotr,double *veold,ITG *nactdof,
	       double *bcont,double *fn));

void FORTRAN(temploadmodal,(double *amta,ITG *namta,ITG *nam,double *ampli,
         double *timemin,double *ttimemin,double *dtime,double *xbounold,
         double *xboun,double *xbounmin,ITG *iamboun,ITG *nboun,
         ITG *nodeboun,ITG *ndirboun,char *amname));
    
void FORTRAN(tiefaccont,(char *lakon,ITG *ipkon,ITG *kon,ITG *ntie,
       char *tieset,ITG *nset,char *set,ITG *istartset,ITG *iendset,
       ITG *ialset,ITG *itiefac,ITG *islavsurf,ITG *islavnode,
       ITG *imastnode,ITG *nslavnode,ITG *nmastnode,ITG *nslavs,
       ITG *nmasts,ITG *ifacecount,ITG *iponoels,ITG *inoels,ITG *ifreenoels,
       ITG *mortar,ITG *ipoface,ITG *nodface,ITG *nk,double *xnoels));   

void tiedcontact(ITG *ntie,char *tieset,ITG *nset,char *set,
               ITG *istartset,ITG *iendset,ITG *ialset,
               char *lakon,ITG *ipkon,ITG *kon,double *tietol,
               ITG *nmpc,ITG *mpcfree,ITG *memmpc_,
               ITG **ipompcp,char **labmpcp,ITG **ikmpcp,ITG **ilmpcp,
               double **fmpcp,ITG **nodempcp,double **coefmpcp,
	       ITG *ithermal,double *co,double *vold,ITG *cfd,
	       ITG *nmpc_,ITG *mi,ITG *nk,ITG *istep,ITG *ikboun,
	       ITG *nboun,char *kind1,char *kind2);
	
void FORTRAN(transformatrix,(double *xab,double *p,double *a));

void FORTRAN(trianeighbor,(ITG *ipe,ITG *ime,ITG *imastop,ITG *ncont,
               ITG *koncont,ITG *ifreeme));

void FORTRAN(triangucont,(ITG *ncont,ITG *ntie,char *tieset,ITG *nset,
          char *set,ITG *istartset,ITG *iendset,ITG *ialset,ITG *itietri,
	  char *lakon,ITG *ipkon,ITG *kon,ITG *koncont,char *kind1,
	  char *kind2,double *co,ITG *nk));

#ifdef BAM
void FORTRAN(uexternaldb,(ITG *lop,ITG *lrestart,double *time,double *dtime,
                          ITG *kstep,ITG *kinc));
#endif

void FORTRAN(ufaceload,(double *co,ITG *ipkon,ITG *kon,char *lakon,
			ITG *nboun,ITG *nodeboun,
                        ITG *nelemload,char *sideload,ITG *nload,
                        ITG *ne,ITG *nk));

void FORTRAN(uinit,());

void FORTRAN(uiter,(ITG *iit));

void FORTRAN(uout,(double *v,ITG *mi,ITG *ithermal,char *filab));

void FORTRAN(updatecont,(ITG *koncont,ITG *ncont,double *co,double *vold,
			 double *cg,double *straight,ITG *mi));

void FORTRAN(updatecontpen,(ITG *koncont,ITG *ncont,double *co,double *vold,
			 double *cg,double *straight,ITG *mi,ITG *imastnode,
                         ITG *nmastnode,double *xmastnor,ITG *ntie,
                         char *tieset,ITG *nset,char *set,ITG *istartset,
                         ITG *iendset,ITG *ialset,ITG *ipkon,char *lakon,
			 ITG *kon,double *cs,ITG *mcs,ITG *ics));

void *u_calloc(size_t num,size_t size,const char *file,const int line,const char *ptr_name);

void *u_free(void* num,const char *file,const int line,const char *ptr_name);

void *u_realloc(void* num,size_t size,const char *file,const int line,const char *ptr_name);

void writeBasisParameter(FILE *f,ITG *istep,ITG *iinc);

void FORTRAN(writeboun,(ITG *nodeboun,ITG *ndirboun,double *xboun,
      char *typeboun,ITG *nboun));

void FORTRAN(writebv,(double *,ITG *));

void FORTRAN(writecvg,(ITG *itep,ITG *iinc,ITG *icutb,ITG *iit,ITG *ne,ITG *ne0,
		       double *ram,double *qam,double *cam,double *uam,
                       ITG *ithermal));

void FORTRAN(writeev,(double *,ITG *,double *,double *));

void FORTRAN(writeevcomplex,(double *eigxx,ITG *nev,double *fmin,double *fmax));

void FORTRAN(writeevcs,(double *,ITG *,ITG *,double *,double *));

void FORTRAN(writeevcscomplex,(double *eigxx,ITG *nev,ITG *nm,double *fmin,
            double *fmax));

void FORTRAN(writehe,(ITG *));

void writeheading(char *jobnamec,char *heading,ITG *nheading);

void FORTRAN(writeim,());

void FORTRAN(writeinput,(char *inpc,ITG *ipoinp,ITG *inp,ITG *nline,ITG *ninp,
                         ITG *ipoinpc));

void FORTRAN(writemac,(double *mac,ITG *nev));

void FORTRAN(writemaccs,(double *mac,ITG *nev,ITG* nm));

void FORTRAN(writempc,(ITG *,ITG *,double *,char *,ITG *));

void FORTRAN(writepf,(double *d,double *bjr,double *bji,double *freq ,
		      ITG *nev,ITG *mode,ITG *nherm));

void FORTRAN(writere,());

void FORTRAN(writesubmatrix,(double *submatrix,ITG *noderetain,
             ITG *ndirretain,ITG *nretain,char *jobnamec));

void FORTRAN(writesummary,(ITG *istep,ITG *j,ITG *icutb,ITG *l,double *ttime,
		   double *time,double *dtime));

void FORTRAN(writesummarydiv,(ITG *istep,ITG *j,ITG *icutb,ITG *l,double *ttime,
		   double *time,double *dtime));

void FORTRAN(writetetmesh,(ITG *kontet,ITG *netet_,double *cotet,
     ITG *nktet,double *field,ITG *nfield));
			
void FORTRAN(writeview,(ITG *ntr,double *adview,double *auview,double *fenv,
            ITG *nzsrad,char *jobnamef));

void FORTRAN(zienzhu,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
		      ITG *ne,double *stn,ITG *ipneigh,ITG *neigh,
		      double *sti,ITG *mi));

void FORTRAN(znaupd,(ITG *ido,char *bmat,ITG *n,char *which,ITG *nev,
	     double *tol,double *resid,ITG *ncv,double *z,ITG *ldz,
	     ITG *iparam,ITG *ipntr,double *workd,double *workl,
	     ITG *lworkl,double *rwork,ITG *info));

void FORTRAN(zneupd,(ITG *rvec,char *howmny,ITG *select,double *d,
	     double *z,ITG *ldz,double *sigma,
             double *workev,char *bmat,ITG *neq,char *which,
	     ITG *nev,double *tol,double *resid,ITG *ncv,double *v,
	     ITG *ldv,ITG *iparam,ITG *ipntr,double *workd,
	     double *workl,ITG *lworkl,double *rwork,ITG *info));

void FORTRAN(getflux,(double *co,ITG *ntmat_,
			double *vold,double *cocon,
			ITG *ncocon,ITG * iset,ITG *istartset,ITG *iendset,
			ITG *ipkon,char *lakon,ITG *kon,ITG *ialset,
			ITG *ielmat,ITG *mi,double * flux));

void FORTRAN(getkdeltatemp,(double *co,ITG *ntmat_,
			double *vold,double *cocon,
			ITG *ncocon,ITG * iset,ITG *istartset,ITG *iendset,
			ITG *ipkon,char *lakon,ITG *kon,ITG *ialset,
			ITG *ielmat,ITG *mi,double * kdelta, double * reftemp));

void FORTRAN(getgausspointscoords,(int setID, double *co,
			ITG * iset,ITG *istartset,ITG *iendset,
			ITG *ipkon,char *lakon,ITG *kon,ITG *ialset));
