/*     CALCULIX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2022 Guido Dhondt                     */

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

/* preCICE header files*/

void nonlingeo_precice(double **co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
                       ITG *ne,
                       ITG *nodeboun, ITG *ndirboun, double *xboun, ITG *nboun,
                       ITG **ipompcp, ITG **nodempcp, double **coefmpcp, char **labmpcp,
                       ITG *nmpc,
                       ITG *nodeforc, ITG *ndirforc, double *xforc, ITG *nforc,
                       ITG **nelemloadp, char **sideloadp, double *xload,
                       ITG *nload, ITG *nactdof,
                       ITG **icolp, ITG *jq, ITG **irowp, ITG *neq, ITG *nzl,
                       ITG *nmethod, ITG **ikmpcp, ITG **ilmpcp, ITG *ikboun,
                       ITG *   ilboun,
                       double *elcon, ITG *nelcon, double *rhcon, ITG *nrhcon,
                       double *alcon, ITG *nalcon, double *alzero, ITG **ielmatp,
                       ITG **ielorienp, ITG *norien, double *orab, ITG *ntmat_,
                       double *t0, double *t1, double *t1old,
                       ITG *ithermal, double *prestr, ITG *iprestr,
                       double **vold, ITG *iperturb, double *sti, ITG *nzs,
                       ITG *kode, char *filab, ITG *idrct,
                       ITG *jmax, ITG *jout, double *timepar,
                       double *eme, double *xbounold,
                       double *xforcold, double *xloadold,
                       double *veold, double *accold,
                       char *amname, double *amta, ITG *namta, ITG *nam,
                       ITG *iamforc, ITG **iamloadp,
                       ITG *iamt1, double *alpha, ITG *iexpl,
                       ITG *iamboun, double *plicon, ITG *nplicon, double *plkcon,
                       ITG *    nplkcon,
                       double **xstatep, ITG *npmat_, ITG *istep, double *ttime,
                       char *matname, double *qaold, ITG *mi,
                       ITG *isolver, ITG *ncmat_, ITG *nstate_, ITG *iumat,
                       double *cs, ITG *mcs, ITG *nkon, double **ener, ITG *mpcinfo,
                       char *  output,
                       double *shcon, ITG *nshcon, double *cocon, ITG *ncocon,
                       double *physcon, ITG *nflow, double *ctrl,
                       char *set, ITG *nset, ITG *istartset,
                       ITG *iendset, ITG *ialset, ITG *nprint, char *prlab,
                       char *prset, ITG *nener, ITG *ikforc, ITG *ilforc, double *trab,
                       ITG *inotr, ITG *ntrans, double **fmpcp, char *cbody,
                       ITG *ibody, double *xbody, ITG *nbody, double *xbodyold,
                       ITG *ielprop, double *prop, ITG *ntie, char *tieset,
                       ITG *itpamp, ITG *iviewfile, char *jobnamec, double *tietol,
                       ITG *nslavs, double *thicke, ITG *ics,
                       ITG *nintpoint, ITG *mortar, ITG *ifacecount, char *typeboun,
                       ITG **islavsurfp, double **pslavsurfp, double **clearinip,
                       ITG *nmat, double *xmodal, ITG *iaxial, ITG *inext, ITG *nprop,
                       ITG *network, char *orname, double *vel, ITG *nef,
                       double *velo, double *veloo, double *energy, ITG *itempuser,
                       ITG *ipobody, ITG *inewton, double *t0g, double *t1g,
                       ITG *ifreebody, char *preciceParticipantName, char *configFilename);

void linstatic_precice(double *co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
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
	       ITG **ielorienp,ITG *norien,double *orab,ITG *ntmat_,
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
	       ITG *inotr,ITG *ntrans,double *fmpc,ITG *ipobody,ITG *ibody,
	       double *xbody,ITG *nbody,double *xbodyold,double *timepar,
	       double *thicke,char *jobnamec,char *tieset,ITG *ntie,
	       ITG *istep,ITG *nmat,ITG *ielprop,double *prop,char *typeboun,
	       ITG *mortar,ITG *mpcinfo,double *tietol,ITG *ics,
	       char *orname,ITG *itempuser,double *t0g,double *t1g,
	       ITG *jmax, char *preciceParticipantName, char *configFilename);
