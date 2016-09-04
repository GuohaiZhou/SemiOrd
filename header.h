
#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED


#include<iostream>
#include <sstream>
#include<fstream>
#include<string>
#include<cmath>
#include<dlib/matrix.h>
//#include<dlib/matrix_la_abstract.h>
#include<dlib/optimization.h>
#include<dlib/rand.h>

#include <scythestat/rng/mersenne.h>
#include <scythestat/distributions.h>
#include <scythestat/ide.h>
#include <scythestat/la.h>
#include <scythestat/matrix.h>
#include <scythestat/rng.h>
#include <scythestat/smath.h>
#include <scythestat/stat.h>
#include <scythestat/optimize.h>
#include <time.h>

using namespace dlib;
using namespace scythe;
std::ofstream DebugOut ("Debug_info.txt", std::ofstream::out);


const double BIGNUM = std::numeric_limits<double>::infinity();
const double SigLevel = 0.05;
// const double pi=std::asin(1)*2;
// dlib::pi
typedef matrix<double,0,1> dbcolvec;
typedef matrix<double> dbmat;
typedef matrix<int> intmat;
// ############# optimization routines equipped with Convergence check
class bobyqa_implementationConvgFail
    {
        typedef long integer;
        typedef double doublereal;

    public:

        template <
            typename funct,
            typename T,
            typename U
            >
        double find_min (
            const funct& f,
            T& x,
            long npt,
            const U& xl_,
            const U& xu_,
            const double rhobeg,
            const double rhoend,
            const long max_f_evals,  int& ConvgFail
        ) const
        {
            const unsigned long n = x.size();
            const unsigned long w_size = (npt+5)*(npt+n)+3*n*(n+5)/2;
            scoped_ptr<doublereal[]> w(new doublereal[w_size]);

            // make these temporary matrices becuse U might be some
            // kind of matrix_exp that doesn't support taking the address
            // of the first element.
            matrix<double,0,1> xl(xl_);
            matrix<double,0,1> xu(xu_);


            return bobyqa_ (f,
                            x.size(),
                            npt,
                            &x(0),
                            &xl(0),
                            &xu(0),
                            rhobeg,
                            rhoend,
                            max_f_evals,
                            w.get(), ConvgFail );
        }

    private:


        template <typename funct>
        doublereal bobyqa_(
            const funct& calfun,
            const integer n,
            const integer npt,
            doublereal *x,
            const doublereal *xl,
            const doublereal *xu,
            const doublereal rhobeg,
            const doublereal rhoend,
            const integer maxfun,
            doublereal *w,  int& ConvgFail
        ) const
        {

            /* System generated locals */
            integer i__1;
            doublereal d__1, d__2;

            /* Local variables */
            integer j, id_, np, iw, igo, ihq, ixb, ixa, ifv, isl, jsl, ipq, ivl, ixn, ixo, ixp, isu, jsu, ndim;
            doublereal temp, zero;
            integer ibmat, izmat;


            /*     This subroutine seeks the least value of a function of many variables, */
            /*     by applying a trust region method that forms quadratic models by */
            /*     interpolation. There is usually some freedom in the interpolation */
            /*     conditions, which is taken up by minimizing the Frobenius norm of */
            /*     the change to the second derivative of the model, beginning with the */
            /*     zero matrix. The values of the variables are constrained by upper and */
            /*     lower bounds. The arguments of the subroutine are as follows. */

            /*     N must be set to the number of variables and must be at least two. */
            /*     NPT is the number of interpolation conditions. Its value must be in */
            /*       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not */
            /*       recommended. */
            /*     Initial values of the variables must be set in X(1),X(2),...,X(N). They */
            /*       will be changed to the values that give the least calculated F. */
            /*     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper */
            /*       bounds, respectively, on X(I). The construction of quadratic models */
            /*       requires XL(I) to be strictly less than XU(I) for each I. Further, */
            /*       the contribution to a model from changes to the I-th variable is */
            /*       damaged severely by rounding errors if XU(I)-XL(I) is too small. */
            /*     RHOBEG and RHOEND must be set to the initial and final values of a trust */
            /*       region radius, so both must be positive with RHOEND no greater than */
            /*       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest */
            /*       expected change to a variable, while RHOEND should indicate the */
            /*       accuracy that is required in the final values of the variables. An */
            /*       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N, */
            /*       is less than 2*RHOBEG. */
            /*     MAXFUN must be set to an upper bound on the number of calls of CALFUN. */
            /*     The array W will be used for working space. Its length must be at least */
            /*       (NPT+5)*(NPT+N)+3*N*(N+5)/2. */

            /* Parameter adjustments */
            --w;
            --xu;
            --xl;
            --x;

            /* Function Body */
            np = n + 1;

            /*     Return if the value of NPT is unacceptable. */
            if (npt < n + 2 || npt > (n + 2) * np / 2) {
    throw bobyqa_failure("Return from BOBYQA because NPT is not in the required interval");
    // ConvgFail=1;  return 0;
                //goto L40;
            }

            /*     Partition the working space array, so that different parts of it can */
            /*     be treated separately during the calculation of BOBYQB. The partition */
            /*     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the */
            /*     space that is taken by the last array in the argument list of BOBYQB. */

            ndim = npt + n;
            ixb = 1;
            ixp = ixb + n;
            ifv = ixp + n * npt;
            ixo = ifv + npt;
            igo = ixo + n;
            ihq = igo + n;
            ipq = ihq + n * np / 2;
            ibmat = ipq + npt;
            izmat = ibmat + ndim * n;
            isl = izmat + npt * (npt - np);
            isu = isl + n;
            ixn = isu + n;
            ixa = ixn + n;
            id_ = ixa + n;
            ivl = id_ + n;
            iw = ivl + ndim;

            /*     Return if there is insufficient space between the bounds. Modify the */
            /*     initial X if necessary in order to avoid conflicts between the bounds */
            /*     and the construction of the first quadratic model. The lower and upper */
            /*     bounds on moves from the updated X are set now, in the ISL and ISU */
            /*     partitions of W, in order to provide useful and exact information about */
            /*     components of X that become within distance RHOBEG from their bounds. */

            zero = 0.;
            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                temp = xu[j] - xl[j];
                if (temp < rhobeg + rhobeg) {
                    throw bobyqa_failure("Return from BOBYQA because one of the differences in x_lower and x_upper is less than 2*rho_begin");
                    //goto L40;
                }
                jsl = isl + j - 1;
                jsu = jsl + n;
                w[jsl] = xl[j] - x[j];
                w[jsu] = xu[j] - x[j];
                if (w[jsl] >= -(rhobeg)) {
                    if (w[jsl] >= zero) {
                        x[j] = xl[j];
                        w[jsl] = zero;
                        w[jsu] = temp;
                    } else {
                        x[j] = xl[j] + rhobeg;
                        w[jsl] = -(rhobeg);
                        /* Computing MAX */
                        d__1 = xu[j] - x[j];
                        w[jsu] = std::max(d__1,rhobeg);
                    }
                } else if (w[jsu] <= rhobeg) {
                    if (w[jsu] <= zero) {
                        x[j] = xu[j];
                        w[jsl] = -temp;
                        w[jsu] = zero;
                    } else {
                        x[j] = xu[j] - rhobeg;
                        /* Computing MIN */
                        d__1 = xl[j] - x[j], d__2 = -(rhobeg);
                        w[jsl] = std::min(d__1,d__2);
                        w[jsu] = rhobeg;
                    }
                }
                /* L30: */
            }

            /*     Make the call of BOBYQB. */

            return bobyqb_(calfun, n, npt, &x[1], &xl[1], &xu[1], rhobeg, rhoend, maxfun, &w[
                    ixb], &w[ixp], &w[ifv], &w[ixo], &w[igo], &w[ihq], &w[ipq], &w[
                    ibmat], &w[izmat], ndim, &w[isl], &w[isu], &w[ixn], &w[ixa], &w[
                    id_], &w[ivl], &w[iw],ConvgFail);
            //L40:
            ;
        } /* bobyqa_ */

    // ----------------------------------------------------------------------------------------

        template <typename funct>
        doublereal bobyqb_(
            const funct& calfun,
            const integer n,
            const integer npt,
            doublereal *x,
            const doublereal *xl,
            const doublereal *xu,
            const doublereal rhobeg,
            const doublereal rhoend,
            const integer maxfun,
            doublereal *xbase,
            doublereal *xpt,
            doublereal *fval,
            doublereal *xopt,
            doublereal *gopt,
            doublereal *hq,
            doublereal *pq,
            doublereal *bmat,
            doublereal *zmat,
            const integer ndim,
            doublereal *sl,
            doublereal *su,
            doublereal *xnew,
            doublereal *xalt,
            doublereal *d__,
            doublereal *vlag,
            doublereal *w,  int& ConvgFail
        ) const
        {
            /* System generated locals */
            integer xpt_dim1, xpt_offset, bmat_dim1, bmat_offset, zmat_dim1,
            zmat_offset, i__1, i__2, i__3;
            doublereal d__1, d__2, d__3, d__4;

            /* Local variables */
            doublereal f = 0;
            integer i__, j, k, ih, nf, jj, nh, ip, jp;
            doublereal dx;
            integer np;
            doublereal den = 0, one = 0, ten = 0, dsq = 0, rho = 0, sum = 0, two = 0, diff = 0, half = 0, beta = 0, gisq = 0;
            integer knew = 0;
            doublereal temp, suma, sumb, bsum, fopt;
            integer kopt = 0, nptm;
            doublereal zero, curv;
            integer ksav;
            doublereal gqsq = 0, dist = 0, sumw = 0, sumz = 0, diffa = 0, diffb = 0, diffc = 0, hdiag = 0;
            integer kbase;
            doublereal alpha = 0, delta = 0, adelt = 0, denom = 0, fsave = 0, bdtol = 0, delsq = 0;
            integer nresc, nfsav;
            doublereal ratio = 0, dnorm = 0, vquad = 0, pqold = 0, tenth = 0;
            integer itest;
            doublereal sumpq, scaden;
            doublereal errbig, cauchy, fracsq, biglsq, densav;
            doublereal bdtest;
            doublereal crvmin, frhosq;
            doublereal distsq;
            integer ntrits;
            doublereal xoptsq;



            /*     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN */
            /*       are identical to the corresponding arguments in SUBROUTINE BOBYQA. */
            /*     XBASE holds a shift of origin that should reduce the contributions */
            /*       from rounding errors to values of the model and Lagrange functions. */
            /*     XPT is a two-dimensional array that holds the coordinates of the */
            /*       interpolation points relative to XBASE. */
            /*     FVAL holds the values of F at the interpolation points. */
            /*     XOPT is set to the displacement from XBASE of the trust region centre. */
            /*     GOPT holds the gradient of the quadratic model at XBASE+XOPT. */
            /*     HQ holds the explicit second derivatives of the quadratic model. */
            /*     PQ contains the parameters of the implicit second derivatives of the */
            /*       quadratic model. */
            /*     BMAT holds the last N columns of H. */
            /*     ZMAT holds the factorization of the leading NPT by NPT submatrix of H, */
            /*       this factorization being ZMAT times ZMAT^T, which provides both the */
            /*       correct rank and positive semi-definiteness. */
            /*     NDIM is the first dimension of BMAT and has the value NPT+N. */
            /*     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively. */
            /*       All the components of every XOPT are going to satisfy the bounds */
            /*       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when */
            /*       XOPT is on a constraint boundary. */
            /*     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the */
            /*       vector of variables for the next call of CALFUN. XNEW also satisfies */
            /*       the SL and SU constraints in the way that has just been mentioned. */
            /*     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW */
            /*       in order to increase the denominator in the updating of UPDATE. */
            /*     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT. */
            /*     VLAG contains the values of the Lagrange functions at a new point X. */
            /*       They are part of a product that requires VLAG to be of length NDIM. */
            /*     W is a one-dimensional array that is used for working space. Its length */
            /*       must be at least 3*NDIM = 3*(NPT+N). */

            /*     Set some constants. */

            /* Parameter adjustments */
            zmat_dim1 = npt;
            zmat_offset = 1 + zmat_dim1;
            zmat -= zmat_offset;
            xpt_dim1 = npt;
            xpt_offset = 1 + xpt_dim1;
            xpt -= xpt_offset;
            --x;
            --xl;
            --xu;
            --xbase;
            --fval;
            --xopt;
            --gopt;
            --hq;
            --pq;
            bmat_dim1 = ndim;
            bmat_offset = 1 + bmat_dim1;
            bmat -= bmat_offset;
            --sl;
            --su;
            --xnew;
            --xalt;
            --d__;
            --vlag;
            --w;

            /* Function Body */
            half = .5;
            one = 1.;
            ten = 10.;
            tenth = .1;
            two = 2.;
            zero = 0.;
            np = n + 1;
            nptm = npt - np;
            nh = n * np / 2;

            /*     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ, */
            /*     BMAT and ZMAT for the first iteration, with the corresponding values of */
            /*     of NF and KOPT, which are the number of calls of CALFUN so far and the */
            /*     index of the interpolation point at the trust region centre. Then the */
            /*     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is */
            /*     less than NPT. GOPT will be updated if KOPT is different from KBASE. */

            prelim_(calfun, n, npt, &x[1], &xl[1], &xu[1], rhobeg, maxfun, &xbase[1],
                    &xpt[xpt_offset], &fval[1], &gopt[1], &hq[1], &pq[1], &bmat[bmat_offset],
                    &zmat[zmat_offset], ndim, &sl[1], &su[1], nf, kopt);
            xoptsq = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                xopt[i__] = xpt[kopt + i__ * xpt_dim1];
                /* L10: */
                /* Computing 2nd power */
                d__1 = xopt[i__];
                xoptsq += d__1 * d__1;
            }
            fsave = fval[1];
            if (nf < npt) {
               // throw bobyqa_failure("Return from BOBYQA because the objective function has been called max_f_evals times.");
ConvgFail=1;  return 0;  //  objective function has been called max_f_evals times.
                //goto L720;
            }
            kbase = 1;

            /*     Complete the settings that are required for the iterative procedure. */

            rho = rhobeg;
            delta = rho;
            nresc = nf;
            ntrits = 0;
            diffa = zero;
            diffb = zero;
            itest = 0;
            nfsav = nf;

            /*     Update GOPT if necessary before the first iteration and after each */
            /*     call of RESCUE that makes a call of CALFUN. */

L20:
            if (kopt != kbase) {
                ih = 0;
                i__1 = n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        ++ih;
                        if (i__ < j) {
                            gopt[j] += hq[ih] * xopt[i__];
                        }
                        /* L30: */
                        gopt[i__] += hq[ih] * xopt[j];
                    }
                }
                if (nf > npt) {
                    i__2 = npt;
                    for (k = 1; k <= i__2; ++k) {
                        temp = zero;
                        i__1 = n;
                        for (j = 1; j <= i__1; ++j) {
                            /* L40: */
                            temp += xpt[k + j * xpt_dim1] * xopt[j];
                        }
                        temp = pq[k] * temp;
                        i__1 = n;
                        for (i__ = 1; i__ <= i__1; ++i__) {
                            /* L50: */
                            gopt[i__] += temp * xpt[k + i__ * xpt_dim1];
                        }
                    }
                }
            }

            /*     Generate the next point in the trust region that provides a small value */
            /*     of the quadratic model subject to the constraints on the variables. */
            /*     The integer NTRITS is set to the number "trust region" iterations that */
            /*     have occurred since the last "alternative" iteration. If the length */
            /*     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to */
            /*     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW. */

L60:
            trsbox_(n, npt, &xpt[xpt_offset], &xopt[1], &gopt[1], &hq[1], &pq[1], &sl[1],
                    &su[1], delta, &xnew[1], &d__[1], &w[1], &w[np], &w[np + n],
                    &w[np + (n << 1)], &w[np + n * 3], &dsq, &crvmin);
            /* Computing MIN */
            d__1 = delta, d__2 = std::sqrt(dsq);
            dnorm = std::min(d__1,d__2);
            if (dnorm < half * rho) {
                ntrits = -1;
                /* Computing 2nd power */
                d__1 = ten * rho;
                distsq = d__1 * d__1;
                if (nf <= nfsav + 2) {
                    goto L650;
                }

                /*     The following choice between labels 650 and 680 depends on whether or */
                /*     not our work with the current RHO seems to be complete. Either RHO is */
                /*     decreased or termination occurs if the errors in the quadratic model at */
                /*     the last three interpolation points compare favourably with predictions */
                /*     of likely improvements to the model within distance HALF*RHO of XOPT. */

                /* Computing MAX */
                d__1 = std::max(diffa,diffb);
                errbig = std::max(d__1,diffc);
                frhosq = rho * .125 * rho;
                if (crvmin > zero && errbig > frhosq * crvmin) {
                    goto L650;
                }
                bdtol = errbig / rho;
                i__1 = n;
                for (j = 1; j <= i__1; ++j) {
                    bdtest = bdtol;
                    if (xnew[j] == sl[j]) {
                        bdtest = w[j];
                    }
                    if (xnew[j] == su[j]) {
                        bdtest = -w[j];
                    }
                    if (bdtest < bdtol) {
                        curv = hq[(j + j * j) / 2];
                        i__2 = npt;
                        for (k = 1; k <= i__2; ++k) {
                            /* L70: */
                            /* Computing 2nd power */
                            d__1 = xpt[k + j * xpt_dim1];
                            curv += pq[k] * (d__1 * d__1);
                        }
                        bdtest += half * curv * rho;
                        if (bdtest < bdtol) {
                            goto L650;
                        }
                    }
                    /* L80: */
                }
                goto L680;
            }
            ++ntrits;

            /*     Severe cancellation is likely to occur if XOPT is too far from XBASE. */
            /*     If the following test holds, then XBASE is shifted so that XOPT becomes */
            /*     zero. The appropriate changes are made to BMAT and to the second */
            /*     derivatives of the current model, beginning with the changes to BMAT */
            /*     that do not depend on ZMAT. VLAG is used temporarily for working space. */

L90:
            if (dsq <= xoptsq * .001) {
                fracsq = xoptsq * .25;
                sumpq = zero;
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    sumpq += pq[k];
                    sum = -half * xoptsq;
                    i__2 = n;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        /* L100: */
                        sum += xpt[k + i__ * xpt_dim1] * xopt[i__];
                    }
                    w[npt + k] = sum;
                    temp = fracsq - half * sum;
                    i__2 = n;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        w[i__] = bmat[k + i__ * bmat_dim1];
                        vlag[i__] = sum * xpt[k + i__ * xpt_dim1] + temp * xopt[i__];
                        ip = npt + i__;
                        i__3 = i__;
                        for (j = 1; j <= i__3; ++j) {
                            /* L110: */
                            bmat[ip + j * bmat_dim1] = bmat[ip + j * bmat_dim1] + w[
                                i__] * vlag[j] + vlag[i__] * w[j];
                        }
                    }
                }

                /*     Then the revisions of BMAT that depend on ZMAT are calculated. */

                i__3 = nptm;
                for (jj = 1; jj <= i__3; ++jj) {
                    sumz = zero;
                    sumw = zero;
                    i__2 = npt;
                    for (k = 1; k <= i__2; ++k) {
                        sumz += zmat[k + jj * zmat_dim1];
                        vlag[k] = w[npt + k] * zmat[k + jj * zmat_dim1];
                        /* L120: */
                        sumw += vlag[k];
                    }
                    i__2 = n;
                    for (j = 1; j <= i__2; ++j) {
                        sum = (fracsq * sumz - half * sumw) * xopt[j];
                        i__1 = npt;
                        for (k = 1; k <= i__1; ++k) {
                            /* L130: */
                            sum += vlag[k] * xpt[k + j * xpt_dim1];
                        }
                        w[j] = sum;
                        i__1 = npt;
                        for (k = 1; k <= i__1; ++k) {
                            /* L140: */
                            bmat[k + j * bmat_dim1] += sum * zmat[k + jj * zmat_dim1];
                        }
                    }
                    i__1 = n;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        ip = i__ + npt;
                        temp = w[i__];
                        i__2 = i__;
                        for (j = 1; j <= i__2; ++j) {
                            /* L150: */
                            bmat[ip + j * bmat_dim1] += temp * w[j];
                        }
                    }
                }

                /*     The following instructions complete the shift, including the changes */
                /*     to the second derivative parameters of the quadratic model. */

                ih = 0;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    w[j] = -half * sumpq * xopt[j];
                    i__1 = npt;
                    for (k = 1; k <= i__1; ++k) {
                        w[j] += pq[k] * xpt[k + j * xpt_dim1];
                        /* L160: */
                        xpt[k + j * xpt_dim1] -= xopt[j];
                    }
                    i__1 = j;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        ++ih;
                        hq[ih] = hq[ih] + w[i__] * xopt[j] + xopt[i__] * w[j];
                        /* L170: */
                        bmat[npt + i__ + j * bmat_dim1] = bmat[npt + j + i__ *
                            bmat_dim1];
                    }
                }
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    xbase[i__] += xopt[i__];
                    xnew[i__] -= xopt[i__];
                    sl[i__] -= xopt[i__];
                    su[i__] -= xopt[i__];
                    /* L180: */
                    xopt[i__] = zero;
                }
                xoptsq = zero;
            }
            if (ntrits == 0) {
                goto L210;
            }
            goto L230;

            /*     XBASE is also moved to XOPT by a call of RESCUE. This calculation is */
            /*     more expensive than the previous shift, because new matrices BMAT and */
            /*     ZMAT are generated from scratch, which may include the replacement of */
            /*     interpolation points whose positions seem to be causing near linear */
            /*     dependence in the interpolation conditions. Therefore RESCUE is called */
            /*     only if rounding errors have reduced by at least a factor of two the */
            /*     denominator of the formula for updating the H matrix. It provides a */
            /*     useful safeguard, but is not invoked in most applications of BOBYQA. */

L190:
            nfsav = nf;
            kbase = kopt;
            rescue_(calfun, n, npt, &xl[1], &xu[1], maxfun, &xbase[1], &xpt[
                    xpt_offset], &fval[1], &xopt[1], &gopt[1], &hq[1], &pq[1], &bmat[
                    bmat_offset], &zmat[zmat_offset], ndim, &sl[1], &su[1], nf, delta,
                    kopt, &vlag[1], &w[1], &w[n + np], &w[ndim + np]);

            /*     XOPT is updated now in case the branch below to label 720 is taken. */
            /*     Any updating of GOPT occurs after the branch below to label 20, which */
            /*     leads to a trust region iteration as does the branch to label 60. */

            xoptsq = zero;
            if (kopt != kbase) {
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    xopt[i__] = xpt[kopt + i__ * xpt_dim1];
                    /* L200: */
                    /* Computing 2nd power */
                    d__1 = xopt[i__];
                    xoptsq += d__1 * d__1;
                }
            }
            if (nf < 0) {
                nf = maxfun;
                // throw bobyqa_failure("Return from BOBYQA because the objective function has been called max_f_evals times.");
     ConvgFail=1;  return 0;
                //goto L720;
            }
            nresc = nf;
            if (nfsav < nf) {
                nfsav = nf;
                goto L20;
            }
            if (ntrits > 0) {
                goto L60;
            }

            /*     Pick two alternative vectors of variables, relative to XBASE, that */
            /*     are suitable as new positions of the KNEW-th interpolation point. */
            /*     Firstly, XNEW is set to the point on a line through XOPT and another */
            /*     interpolation point that minimizes the predicted value of the next */
            /*     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL */
            /*     and SU bounds. Secondly, XALT is set to the best feasible point on */
            /*     a constrained version of the Cauchy step of the KNEW-th Lagrange */
            /*     function, the corresponding value of the square of this function */
            /*     being returned in CAUCHY. The choice between these alternatives is */
            /*     going to be made when the denominator is calculated. */

L210:
            altmov_(n, npt, &xpt[xpt_offset], &xopt[1], &bmat[bmat_offset], &zmat[zmat_offset],
                    ndim, &sl[1], &su[1], kopt, knew, adelt, &xnew[1],
                    &xalt[1], alpha, cauchy, &w[1], &w[np], &w[ndim + 1]);
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                /* L220: */
                d__[i__] = xnew[i__] - xopt[i__];
            }

            /*     Calculate VLAG and BETA for the current choice of D. The scalar */
            /*     product of D with XPT(K,.) is going to be held in W(NPT+K) for */
            /*     use when VQUAD is calculated. */

L230:
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                suma = zero;
                sumb = zero;
                sum = zero;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    suma += xpt[k + j * xpt_dim1] * d__[j];
                    sumb += xpt[k + j * xpt_dim1] * xopt[j];
                    /* L240: */
                    sum += bmat[k + j * bmat_dim1] * d__[j];
                }
                w[k] = suma * (half * suma + sumb);
                vlag[k] = sum;
                /* L250: */
                w[npt + k] = suma;
            }
            beta = zero;
            i__1 = nptm;
            for (jj = 1; jj <= i__1; ++jj) {
                sum = zero;
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    /* L260: */
                    sum += zmat[k + jj * zmat_dim1] * w[k];
                }
                beta -= sum * sum;
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    /* L270: */
                    vlag[k] += sum * zmat[k + jj * zmat_dim1];
                }
            }
            dsq = zero;
            bsum = zero;
            dx = zero;
            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                /* Computing 2nd power */
                d__1 = d__[j];
                dsq += d__1 * d__1;
                sum = zero;
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    /* L280: */
                    sum += w[k] * bmat[k + j * bmat_dim1];
                }
                bsum += sum * d__[j];
                jp = npt + j;
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    /* L290: */
                    sum += bmat[jp + i__ * bmat_dim1] * d__[i__];
                }
                vlag[jp] = sum;
                bsum += sum * d__[j];
                /* L300: */
                dx += d__[j] * xopt[j];
            }
            beta = dx * dx + dsq * (xoptsq + dx + dx + half * dsq) + beta - bsum;
            vlag[kopt] += one;

            /*     If NTRITS is zero, the denominator may be increased by replacing */
            /*     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if */
            /*     rounding errors have damaged the chosen denominator. */

            if (ntrits == 0) {
                /* Computing 2nd power */
                d__1 = vlag[knew];
                denom = d__1 * d__1 + alpha * beta;
                if (denom < cauchy && cauchy > zero) {
                    i__2 = n;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        xnew[i__] = xalt[i__];
                        /* L310: */
                        d__[i__] = xnew[i__] - xopt[i__];
                    }
                    cauchy = zero;
                    goto L230;
                }
                /* Computing 2nd power */
                d__1 = vlag[knew];
                if (denom <= half * (d__1 * d__1)) {
                    if (nf > nresc) {
                        goto L190;
                    }
          //  throw bobyqa_failure("Return from BOBYQA because of much cancellation in a denominator.");
          ConvgFail=1;  return 0;
                    //goto L720;
                }

                /*     Alternatively, if NTRITS is positive, then set KNEW to the index of */
                /*     the next interpolation point to be deleted to make room for a trust */
                /*     region step. Again RESCUE may be called if rounding errors have damaged */
                /*     the chosen denominator, which is the reason for attempting to select */
                /*     KNEW before calculating the next value of the objective function. */

            } else {
                delsq = delta * delta;
                scaden = zero;
                biglsq = zero;
                knew = 0;
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    if (k == kopt) {
                        goto L350;
                    }
                    hdiag = zero;
                    i__1 = nptm;
                    for (jj = 1; jj <= i__1; ++jj) {
                        /* L330: */
                        /* Computing 2nd power */
                        d__1 = zmat[k + jj * zmat_dim1];
                        hdiag += d__1 * d__1;
                    }
                    /* Computing 2nd power */
                    d__1 = vlag[k];
                    den = beta * hdiag + d__1 * d__1;
                    distsq = zero;
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j) {
                        /* L340: */
                        /* Computing 2nd power */
                        d__1 = xpt[k + j * xpt_dim1] - xopt[j];
                        distsq += d__1 * d__1;
                    }
                    /* Computing MAX */
                    /* Computing 2nd power */
                    d__3 = distsq / delsq;
                    d__1 = one, d__2 = d__3 * d__3;
                    temp = std::max(d__1,d__2);
                    if (temp * den > scaden) {
                        scaden = temp * den;
                        knew = k;
                        denom = den;
                    }
                    /* Computing MAX */
                    /* Computing 2nd power */
                    d__3 = vlag[k];
                    d__1 = biglsq, d__2 = temp * (d__3 * d__3);
                    biglsq = std::max(d__1,d__2);
L350:
                    ;
                }
                if (scaden <= half * biglsq) {
                    if (nf > nresc) {
                        goto L190;
                    }
  //  throw bobyqa_failure("Return from BOBYQA because of much cancellation in a denominator.");
                    //goto L720;
   ConvgFail=1;  return 0;
                }
            }

            /*     Put the variables for the next calculation of the objective function */
            /*       in XNEW, with any adjustments for the bounds. */


            /*     Calculate the value of the objective function at XBASE+XNEW, unless */
            /*       the limit on the number of calculations of F has been reached. */

L360:
            i__2 = n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                /* Computing MIN */
                /* Computing MAX */
                d__3 = xl[i__], d__4 = xbase[i__] + xnew[i__];
                d__1 = std::max(d__3,d__4), d__2 = xu[i__];
                x[i__] = std::min(d__1,d__2);
                if (xnew[i__] == sl[i__]) {
                    x[i__] = xl[i__];
                }
                if (xnew[i__] == su[i__]) {
                    x[i__] = xu[i__];
                }
                /* L380: */
            }
            if (nf >= maxfun) {
              //  throw bobyqa_failure("Return from BOBYQA because the objective function has been called max_f_evals times.");
	    ConvgFail=1;  return 0;
                //goto L720;
            }
            ++nf;
            f = calfun(mat(&x[1], n));
            if (ntrits == -1) {
                fsave = f;
                goto L720;
            }

            /*     Use the quadratic model to predict the change in F due to the step D, */
            /*       and set DIFF to the error of this prediction. */

            fopt = fval[kopt];
            vquad = zero;
            ih = 0;
            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                vquad += d__[j] * gopt[j];
                i__1 = j;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    ++ih;
                    temp = d__[i__] * d__[j];
                    if (i__ == j) {
                        temp = half * temp;
                    }
                    /* L410: */
                    vquad += hq[ih] * temp;
                }
            }
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                /* L420: */
                /* Computing 2nd power */
                d__1 = w[npt + k];
                vquad += half * pq[k] * (d__1 * d__1);
            }
            diff = f - fopt - vquad;
            diffc = diffb;
            diffb = diffa;
            diffa = std::abs(diff);
            if (dnorm > rho) {
                nfsav = nf;
            }

            /*     Pick the next value of DELTA after a trust region step. */

            if (ntrits > 0) {
                if (vquad >= zero) {
  //  throw ("Return from BOBYQA because a trust region step has failed to reduce Q.");
                    //goto L720;
 { return 0; ConvgFail=1; }
                }
                ratio = (f - fopt) / vquad;
                if (ratio <= tenth) {
                    /* Computing MIN */
                    d__1 = half * delta;
                    delta = std::min(d__1,dnorm);
                } else if (ratio <= .7) {
                    /* Computing MAX */
                    d__1 = half * delta;
                    delta = std::max(d__1,dnorm);
                } else {
                    /* Computing MAX */
                    d__1 = half * delta, d__2 = dnorm + dnorm;
                    delta = std::max(d__1,d__2);
                }
                if (delta <= rho * 1.5) {
                    delta = rho;
                }

                /*     Recalculate KNEW and DENOM if the new F is less than FOPT. */

                if (f < fopt) {
                    ksav = knew;
                    densav = denom;
                    delsq = delta * delta;
                    scaden = zero;
                    biglsq = zero;
                    knew = 0;
                    i__1 = npt;
                    for (k = 1; k <= i__1; ++k) {
                        hdiag = zero;
                        i__2 = nptm;
                        for (jj = 1; jj <= i__2; ++jj) {
                            /* L440: */
                            /* Computing 2nd power */
                            d__1 = zmat[k + jj * zmat_dim1];
                            hdiag += d__1 * d__1;
                        }
                        /* Computing 2nd power */
                        d__1 = vlag[k];
                        den = beta * hdiag + d__1 * d__1;
                        distsq = zero;
                        i__2 = n;
                        for (j = 1; j <= i__2; ++j) {
                            /* L450: */
                            /* Computing 2nd power */
                            d__1 = xpt[k + j * xpt_dim1] - xnew[j];
                            distsq += d__1 * d__1;
                        }
                        /* Computing MAX */
                        /* Computing 2nd power */
                        d__3 = distsq / delsq;
                        d__1 = one, d__2 = d__3 * d__3;
                        temp = std::max(d__1,d__2);
                        if (temp * den > scaden) {
                            scaden = temp * den;
                            knew = k;
                            denom = den;
                        }
                        /* L460: */
                        /* Computing MAX */
                        /* Computing 2nd power */
                        d__3 = vlag[k];
                        d__1 = biglsq, d__2 = temp * (d__3 * d__3);
                        biglsq = std::max(d__1,d__2);
                    }
                    if (scaden <= half * biglsq) {
                        knew = ksav;
                        denom = densav;
                    }
                }
            }

            /*     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be */
            /*     moved. Also update the second derivative terms of the model. */

            update_(n, npt, &bmat[bmat_offset], &zmat[zmat_offset], ndim, &vlag[1],
                    beta, denom, knew, &w[1]);
            ih = 0;
            pqold = pq[knew];
            pq[knew] = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                temp = pqold * xpt[knew + i__ * xpt_dim1];
                i__2 = i__;
                for (j = 1; j <= i__2; ++j) {
                    ++ih;
                    /* L470: */
                    hq[ih] += temp * xpt[knew + j * xpt_dim1];
                }
            }
            i__2 = nptm;
            for (jj = 1; jj <= i__2; ++jj) {
                temp = diff * zmat[knew + jj * zmat_dim1];
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    /* L480: */
                    pq[k] += temp * zmat[k + jj * zmat_dim1];
                }
            }

            /*     Include the new interpolation point, and make the changes to GOPT at */
            /*     the old XOPT that are caused by the updating of the quadratic model. */

            fval[knew] = f;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                xpt[knew + i__ * xpt_dim1] = xnew[i__];
                /* L490: */
                w[i__] = bmat[knew + i__ * bmat_dim1];
            }
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                suma = zero;
                i__2 = nptm;
                for (jj = 1; jj <= i__2; ++jj) {
                    /* L500: */
                    suma += zmat[knew + jj * zmat_dim1] * zmat[k + jj * zmat_dim1];
                }
                sumb = zero;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    /* L510: */
                    sumb += xpt[k + j * xpt_dim1] * xopt[j];
                }
                temp = suma * sumb;
                i__2 = n;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    /* L520: */
                    w[i__] += temp * xpt[k + i__ * xpt_dim1];
                }
            }
            i__2 = n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                /* L530: */
                gopt[i__] += diff * w[i__];
            }

            /*     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT. */

            if (f < fopt) {
                kopt = knew;
                xoptsq = zero;
                ih = 0;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    xopt[j] = xnew[j];
                    /* Computing 2nd power */
                    d__1 = xopt[j];
                    xoptsq += d__1 * d__1;
                    i__1 = j;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        ++ih;
                        if (i__ < j) {
                            gopt[j] += hq[ih] * d__[i__];
                        }
                        /* L540: */
                        gopt[i__] += hq[ih] * d__[j];
                    }
                }
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    temp = zero;
                    i__2 = n;
                    for (j = 1; j <= i__2; ++j) {
                        /* L550: */
                        temp += xpt[k + j * xpt_dim1] * d__[j];
                    }
                    temp = pq[k] * temp;
                    i__2 = n;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        /* L560: */
                        gopt[i__] += temp * xpt[k + i__ * xpt_dim1];
                    }
                }
            }

            /*     Calculate the parameters of the least Frobenius norm interpolant to */
            /*     the current data, the gradient of this interpolant at XOPT being put */
            /*     into VLAG(NPT+I), I=1,2,...,N. */

            if (ntrits > 0) {
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    vlag[k] = fval[k] - fval[kopt];
                    /* L570: */
                    w[k] = zero;
                }
                i__2 = nptm;
                for (j = 1; j <= i__2; ++j) {
                    sum = zero;
                    i__1 = npt;
                    for (k = 1; k <= i__1; ++k) {
                        /* L580: */
                        sum += zmat[k + j * zmat_dim1] * vlag[k];
                    }
                    i__1 = npt;
                    for (k = 1; k <= i__1; ++k) {
                        /* L590: */
                        w[k] += sum * zmat[k + j * zmat_dim1];
                    }
                }
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    sum = zero;
                    i__2 = n;
                    for (j = 1; j <= i__2; ++j) {
                        /* L600: */
                        sum += xpt[k + j * xpt_dim1] * xopt[j];
                    }
                    w[k + npt] = w[k];
                    /* L610: */
                    w[k] = sum * w[k];
                }
                gqsq = zero;
                gisq = zero;
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    sum = zero;
                    i__2 = npt;
                    for (k = 1; k <= i__2; ++k) {
                        /* L620: */
                        sum = sum + bmat[k + i__ * bmat_dim1] * vlag[k] + xpt[k + i__
                            * xpt_dim1] * w[k];
                    }
                    if (xopt[i__] == sl[i__]) {
                        /* Computing MIN */
                        d__2 = zero, d__3 = gopt[i__];
                        /* Computing 2nd power */
                        d__1 = std::min(d__2,d__3);
                        gqsq += d__1 * d__1;
                        /* Computing 2nd power */
                        d__1 = std::min(zero,sum);
                        gisq += d__1 * d__1;
                    } else if (xopt[i__] == su[i__]) {
                        /* Computing MAX */
                        d__2 = zero, d__3 = gopt[i__];
                        /* Computing 2nd power */
                        d__1 = std::max(d__2,d__3);
                        gqsq += d__1 * d__1;
                        /* Computing 2nd power */
                        d__1 = std::max(zero,sum);
                        gisq += d__1 * d__1;
                    } else {
                        /* Computing 2nd power */
                        d__1 = gopt[i__];
                        gqsq += d__1 * d__1;
                        gisq += sum * sum;
                    }
                    /* L630: */
                    vlag[npt + i__] = sum;
                }

                /*     Test whether to replace the new quadratic model by the least Frobenius */
                /*     norm interpolant, making the replacement if the test is satisfied. */

                ++itest;
                if (gqsq < ten * gisq) {
                    itest = 0;
                }
                if (itest >= 3) {
                    i__1 = std::max(npt,nh);
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        if (i__ <= n) {
                            gopt[i__] = vlag[npt + i__];
                        }
                        if (i__ <= npt) {
                            pq[i__] = w[npt + i__];
                        }
                        if (i__ <= nh) {
                            hq[i__] = zero;
                        }
                        itest = 0;
                        /* L640: */
                    }
                }
            }

            /*     If a trust region step has provided a sufficient decrease in F, then */
            /*     branch for another trust region calculation. The case NTRITS=0 occurs */
            /*     when the new interpolation point was reached by an alternative step. */

            if (ntrits == 0) {
                goto L60;
            }
            if (f <= fopt + tenth * vquad) {
                goto L60;
            }

            /*     Alternatively, find out if the interpolation points are close enough */
            /*       to the best point so far. */

            /* Computing MAX */
            /* Computing 2nd power */
            d__3 = two * delta;
            /* Computing 2nd power */
            d__4 = ten * rho;
            d__1 = d__3 * d__3, d__2 = d__4 * d__4;
            distsq = std::max(d__1,d__2);
L650:
            knew = 0;
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                sum = zero;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    /* L660: */
                    /* Computing 2nd power */
                    d__1 = xpt[k + j * xpt_dim1] - xopt[j];
                    sum += d__1 * d__1;
                }
                if (sum > distsq) {
                    knew = k;
                    distsq = sum;
                }
                /* L670: */
            }

            /*     If KNEW is positive, then ALTMOV finds alternative new positions for */
            /*     the KNEW-th interpolation point within distance ADELT of XOPT. It is */
            /*     reached via label 90. Otherwise, there is a branch to label 60 for */
            /*     another trust region iteration, unless the calculations with the */
            /*     current RHO are complete. */

            if (knew > 0) {
                dist = std::sqrt(distsq);
                if (ntrits == -1) {
                    /* Computing MIN */
                    d__1 = tenth * delta, d__2 = half * dist;
                    delta = std::min(d__1,d__2);
                    if (delta <= rho * 1.5) {
                        delta = rho;
                    }
                }
                ntrits = 0;
                /* Computing MAX */
                /* Computing MIN */
                d__2 = tenth * dist;
                d__1 = std::min(d__2,delta);
                adelt = std::max(d__1,rho);
                dsq = adelt * adelt;
                goto L90;
            }
            if (ntrits == -1) {
                goto L680;
            }
            if (ratio > zero) {
                goto L60;
            }
            if (std::max(delta,dnorm) > rho) {
                goto L60;
            }

            /*     The calculations with the current value of RHO are complete. Pick the */
            /*       next values of RHO and DELTA. */

L680:
            if (rho > rhoend) {
                delta = half * rho;
                ratio = rho / rhoend;
                if (ratio <= 16.) {
                    rho = rhoend;
                } else if (ratio <= 250.) {
                    rho = std::sqrt(ratio) * rhoend;
                } else {
                    rho = tenth * rho;
                }
                delta = std::max(delta,rho);
                ntrits = 0;
                nfsav = nf;
                goto L60;
            }

            /*     Return from the calculation, after another Newton-Raphson step, if */
            /*       it is too short to have been tried before. */

            if (ntrits == -1) {
                goto L360;
            }
L720:
            if (fval[kopt] <= fsave) {
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    /* Computing MIN */
                    /* Computing MAX */
                    d__3 = xl[i__], d__4 = xbase[i__] + xopt[i__];
                    d__1 = std::max(d__3,d__4), d__2 = xu[i__];
                    x[i__] = std::min(d__1,d__2);
                    if (xopt[i__] == sl[i__]) {
                        x[i__] = xl[i__];
                    }
                    if (xopt[i__] == su[i__]) {
                        x[i__] = xu[i__];
                    }
                    /* L730: */
                }
                f = fval[kopt];
            }

            return f;
        } /* bobyqb_ */

    // ----------------------------------------------------------------------------------------

        void altmov_(
            const integer n,
            const integer npt,
            const doublereal *xpt,
            const doublereal *xopt,
            const doublereal *bmat,
            const doublereal *zmat,
            const integer ndim,
            const doublereal *sl,
            const doublereal *su,
            const integer kopt,
            const integer knew,
            const doublereal adelt,
            doublereal *xnew,
            doublereal *xalt,
            doublereal& alpha,
            doublereal& cauchy,
            doublereal *glag,
            doublereal *hcol,
            doublereal *w
        ) const
        {
            /* System generated locals */
            integer xpt_dim1, xpt_offset, bmat_dim1, bmat_offset, zmat_dim1,
            zmat_offset, i__1, i__2;
            doublereal d__1, d__2, d__3, d__4;


            /* Local variables */
            integer i__, j, k;
            doublereal ha, gw, one, diff, half;
            integer ilbd, isbd;
            doublereal slbd;
            integer iubd;
            doublereal vlag, subd, temp;
            integer ksav = 0;
            doublereal step = 0, zero = 0, curv = 0;
            integer iflag;
            doublereal scale = 0, csave = 0, tempa = 0, tempb = 0, tempd = 0, const__ = 0, sumin = 0,
                       ggfree = 0;
            integer ibdsav = 0;
            doublereal dderiv = 0, bigstp = 0, predsq = 0, presav = 0, distsq = 0, stpsav = 0, wfixsq = 0, wsqsav = 0;


            /*     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have */
            /*       the same meanings as the corresponding arguments of BOBYQB. */
            /*     KOPT is the index of the optimal interpolation point. */
            /*     KNEW is the index of the interpolation point that is going to be moved. */
            /*     ADELT is the current trust region bound. */
            /*     XNEW will be set to a suitable new position for the interpolation point */
            /*       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region */
            /*       bounds and it should provide a large denominator in the next call of */
            /*       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the */
            /*       straight lines through XOPT and another interpolation point. */
            /*     XALT also provides a large value of the modulus of the KNEW-th Lagrange */
            /*       function subject to the constraints that have been mentioned, its main */
            /*       difference from XNEW being that XALT-XOPT is a constrained version of */
            /*       the Cauchy step within the trust region. An exception is that XALT is */
            /*       not calculated if all components of GLAG (see below) are zero. */
            /*     ALPHA will be set to the KNEW-th diagonal element of the H matrix. */
            /*     CAUCHY will be set to the square of the KNEW-th Lagrange function at */
            /*       the step XALT-XOPT from XOPT for the vector XALT that is returned, */
            /*       except that CAUCHY is set to zero if XALT is not calculated. */
            /*     GLAG is a working space vector of length N for the gradient of the */
            /*       KNEW-th Lagrange function at XOPT. */
            /*     HCOL is a working space vector of length NPT for the second derivative */
            /*       coefficients of the KNEW-th Lagrange function. */
            /*     W is a working space vector of length 2N that is going to hold the */
            /*       constrained Cauchy step from XOPT of the Lagrange function, followed */
            /*       by the downhill version of XALT when the uphill step is calculated. */

            /*     Set the first NPT components of W to the leading elements of the */
            /*     KNEW-th column of the H matrix. */

            /* Parameter adjustments */
            zmat_dim1 = npt;
            zmat_offset = 1 + zmat_dim1;
            zmat -= zmat_offset;
            xpt_dim1 = npt;
            xpt_offset = 1 + xpt_dim1;
            xpt -= xpt_offset;
            --xopt;
            bmat_dim1 = ndim;
            bmat_offset = 1 + bmat_dim1;
            bmat -= bmat_offset;
            --sl;
            --su;
            --xnew;
            --xalt;
            --glag;
            --hcol;
            --w;

            /* Function Body */
            half = .5;
            one = 1.;
            zero = 0.;
            const__ = one + std::sqrt(2.);
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                /* L10: */
                hcol[k] = zero;
            }
            i__1 = npt - n - 1;
            for (j = 1; j <= i__1; ++j) {
                temp = zmat[knew + j * zmat_dim1];
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    /* L20: */
                    hcol[k] += temp * zmat[k + j * zmat_dim1];
                }
            }
            alpha = hcol[knew];
            ha = half * alpha;

            /*     Calculate the gradient of the KNEW-th Lagrange function at XOPT. */

            i__2 = n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                /* L30: */
                glag[i__] = bmat[knew + i__ * bmat_dim1];
            }
            i__2 = npt;
            for (k = 1; k <= i__2; ++k) {
                temp = zero;
                i__1 = n;
                for (j = 1; j <= i__1; ++j) {
                    /* L40: */
                    temp += xpt[k + j * xpt_dim1] * xopt[j];
                }
                temp = hcol[k] * temp;
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    /* L50: */
                    glag[i__] += temp * xpt[k + i__ * xpt_dim1];
                }
            }

            /*     Search for a large denominator along the straight lines through XOPT */
            /*     and another interpolation point. SLBD and SUBD will be lower and upper */
            /*     bounds on the step along each of these lines in turn. PREDSQ will be */
            /*     set to the square of the predicted denominator for each line. PRESAV */
            /*     will be set to the largest admissible value of PREDSQ that occurs. */

            presav = zero;
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                if (k == kopt) {
                    goto L80;
                }
                dderiv = zero;
                distsq = zero;
                i__2 = n;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp = xpt[k + i__ * xpt_dim1] - xopt[i__];
                    dderiv += glag[i__] * temp;
                    /* L60: */
                    distsq += temp * temp;
                }
                subd = adelt / std::sqrt(distsq);
                slbd = -subd;
                ilbd = 0;
                iubd = 0;
                sumin = std::min(one,subd);

                /*     Revise SLBD and SUBD if necessary because of the bounds in SL and SU. */

                i__2 = n;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp = xpt[k + i__ * xpt_dim1] - xopt[i__];
                    if (temp > zero) {
                        if (slbd * temp < sl[i__] - xopt[i__]) {
                            slbd = (sl[i__] - xopt[i__]) / temp;
                            ilbd = -i__;
                        }
                        if (subd * temp > su[i__] - xopt[i__]) {
                            /* Computing MAX */
                            d__1 = sumin, d__2 = (su[i__] - xopt[i__]) / temp;
                            subd = std::max(d__1,d__2);
                            iubd = i__;
                        }
                    } else if (temp < zero) {
                        if (slbd * temp > su[i__] - xopt[i__]) {
                            slbd = (su[i__] - xopt[i__]) / temp;
                            ilbd = i__;
                        }
                        if (subd * temp < sl[i__] - xopt[i__]) {
                            /* Computing MAX */
                            d__1 = sumin, d__2 = (sl[i__] - xopt[i__]) / temp;
                            subd = std::max(d__1,d__2);
                            iubd = -i__;
                        }
                    }
                    /* L70: */
                }

                /*     Seek a large modulus of the KNEW-th Lagrange function when the index */
                /*     of the other interpolation point on the line through XOPT is KNEW. */

                if (k == knew) {
                    diff = dderiv - one;
                    step = slbd;
                    vlag = slbd * (dderiv - slbd * diff);
                    isbd = ilbd;
                    temp = subd * (dderiv - subd * diff);
                    if (std::abs(temp) > std::abs(vlag)) {
                        step = subd;
                        vlag = temp;
                        isbd = iubd;
                    }
                    tempd = half * dderiv;
                    tempa = tempd - diff * slbd;
                    tempb = tempd - diff * subd;
                    if (tempa * tempb < zero) {
                        temp = tempd * tempd / diff;
                        if (std::abs(temp) > std::abs(vlag)) {
                            step = tempd / diff;
                            vlag = temp;
                            isbd = 0;
                        }
                    }

                    /*     Search along each of the other lines through XOPT and another point. */

                } else {
                    step = slbd;
                    vlag = slbd * (one - slbd);
                    isbd = ilbd;
                    temp = subd * (one - subd);
                    if (std::abs(temp) > std::abs(vlag)) {
                        step = subd;
                        vlag = temp;
                        isbd = iubd;
                    }
                    if (subd > half) {
                        if (std::abs(vlag) < .25) {
                            step = half;
                            vlag = .25;
                            isbd = 0;
                        }
                    }
                    vlag *= dderiv;
                }

                /*     Calculate PREDSQ for the current line search and maintain PRESAV. */

                temp = step * (one - step) * distsq;
                predsq = vlag * vlag * (vlag * vlag + ha * temp * temp);
                if (predsq > presav) {
                    presav = predsq;
                    ksav = k;
                    stpsav = step;
                    ibdsav = isbd;
                }
L80:
                ;
            }

            /*     Construct XNEW in a way that satisfies the bound constraints exactly. */

            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                temp = xopt[i__] + stpsav * (xpt[ksav + i__ * xpt_dim1] - xopt[i__]);
                /* L90: */
                /* Computing MAX */
                /* Computing MIN */
                d__3 = su[i__];
                d__1 = sl[i__], d__2 = std::min(d__3,temp);
                xnew[i__] = std::max(d__1,d__2);
            }
            if (ibdsav < 0) {
                xnew[-ibdsav] = sl[-ibdsav];
            }
            if (ibdsav > 0) {
                xnew[ibdsav] = su[ibdsav];
            }

            /*     Prepare for the iterative method that assembles the constrained Cauchy */
            /*     step in W. The sum of squares of the fixed components of W is formed in */
            /*     WFIXSQ, and the free components of W are set to BIGSTP. */

            bigstp = adelt + adelt;
            iflag = 0;
L100:
            wfixsq = zero;
            ggfree = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                w[i__] = zero;
                /* Computing MIN */
                d__1 = xopt[i__] - sl[i__], d__2 = glag[i__];
                tempa = std::min(d__1,d__2);
                /* Computing MAX */
                d__1 = xopt[i__] - su[i__], d__2 = glag[i__];
                tempb = std::max(d__1,d__2);
                if (tempa > zero || tempb < zero) {
                    w[i__] = bigstp;
                    /* Computing 2nd power */
                    d__1 = glag[i__];
                    ggfree += d__1 * d__1;
                }
                /* L110: */
            }
            if (ggfree == zero) {
                cauchy = zero;
                goto L200;
            }

            /*     Investigate whether more components of W can be fixed. */

L120:
            temp = adelt * adelt - wfixsq;
            if (temp > zero) {
                wsqsav = wfixsq;
                step = std::sqrt(temp / ggfree);
                ggfree = zero;
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    if (w[i__] == bigstp) {
                        temp = xopt[i__] - step * glag[i__];
                        if (temp <= sl[i__]) {
                            w[i__] = sl[i__] - xopt[i__];
                            /* Computing 2nd power */
                            d__1 = w[i__];
                            wfixsq += d__1 * d__1;
                        } else if (temp >= su[i__]) {
                            w[i__] = su[i__] - xopt[i__];
                            /* Computing 2nd power */
                            d__1 = w[i__];
                            wfixsq += d__1 * d__1;
                        } else {
                            /* Computing 2nd power */
                            d__1 = glag[i__];
                            ggfree += d__1 * d__1;
                        }
                    }
                    /* L130: */
                }
                if (wfixsq > wsqsav && ggfree > zero) {
                    goto L120;
                }
            }

            /*     Set the remaining free components of W and all components of XALT, */
            /*     except that W may be scaled later. */

            gw = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (w[i__] == bigstp) {
                    w[i__] = -step * glag[i__];
                    /* Computing MAX */
                    /* Computing MIN */
                    d__3 = su[i__], d__4 = xopt[i__] + w[i__];
                    d__1 = sl[i__], d__2 = std::min(d__3,d__4);
                    xalt[i__] = std::max(d__1,d__2);
                } else if (w[i__] == zero) {
                    xalt[i__] = xopt[i__];
                } else if (glag[i__] > zero) {
                    xalt[i__] = sl[i__];
                } else {
                    xalt[i__] = su[i__];
                }
                /* L140: */
                gw += glag[i__] * w[i__];
            }

            /*     Set CURV to the curvature of the KNEW-th Lagrange function along W. */
            /*     Scale W by a factor less than one if that can reduce the modulus of */
            /*     the Lagrange function at XOPT+W. Set CAUCHY to the final value of */
            /*     the square of this function. */

            curv = zero;
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                temp = zero;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    /* L150: */
                    temp += xpt[k + j * xpt_dim1] * w[j];
                }
                /* L160: */
                curv += hcol[k] * temp * temp;
            }
            if (iflag == 1) {
                curv = -curv;
            }
            if (curv > -gw && curv < -const__ * gw) {
                scale = -gw / curv;
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    temp = xopt[i__] + scale * w[i__];
                    /* L170: */
                    /* Computing MAX */
                    /* Computing MIN */
                    d__3 = su[i__];
                    d__1 = sl[i__], d__2 = std::min(d__3,temp);
                    xalt[i__] = std::max(d__1,d__2);
                }
                /* Computing 2nd power */
                d__1 = half * gw * scale;
                cauchy = d__1 * d__1;
            } else {
                /* Computing 2nd power */
                d__1 = gw + half * curv;
                cauchy = d__1 * d__1;
            }

            /*     If IFLAG is zero, then XALT is calculated as before after reversing */
            /*     the sign of GLAG. Thus two XALT vectors become available. The one that */
            /*     is chosen is the one that gives the larger value of CAUCHY. */

            if (iflag == 0) {
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    glag[i__] = -glag[i__];
                    /* L180: */
                    w[n + i__] = xalt[i__];
                }
                csave = cauchy;
                iflag = 1;
                goto L100;
            }
            if (csave > cauchy) {
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    /* L190: */
                    xalt[i__] = w[n + i__];
                }
                cauchy = csave;
            }
L200:
            ;
        } /* altmov_ */

    // ----------------------------------------------------------------------------------------

        template <typename funct>
        void prelim_(
            const funct& calfun,
            const integer n,
            const integer npt,
            doublereal *x,
            const doublereal *xl,
            const doublereal *xu,
            const doublereal rhobeg,
            const integer maxfun,
            doublereal *xbase,
            doublereal *xpt,
            doublereal *fval,
            doublereal *gopt,
            doublereal *hq,
            doublereal *pq,
            doublereal *bmat,
            doublereal *zmat,
            const integer ndim,
            const doublereal *sl,
            const doublereal *su,
            integer& nf,
            integer& kopt
        ) const
        {
            /* System generated locals */
            integer xpt_dim1, xpt_offset, bmat_dim1, bmat_offset, zmat_dim1,
            zmat_offset, i__1, i__2;
            doublereal d__1, d__2, d__3, d__4;


            /* Local variables */
            doublereal f;
            integer i__, j, k, ih, np, nfm;
            doublereal one;
            integer nfx = 0, ipt = 0, jpt = 0;
            doublereal two = 0, fbeg = 0, diff = 0, half = 0, temp = 0, zero = 0, recip = 0, stepa = 0, stepb = 0;
            integer itemp;
            doublereal rhosq;



            /*     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the */
            /*       same as the corresponding arguments in SUBROUTINE BOBYQA. */
            /*     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU */
            /*       are the same as the corresponding arguments in BOBYQB, the elements */
            /*       of SL and SU being set in BOBYQA. */
            /*     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but */
            /*       it is set by PRELIM to the gradient of the quadratic model at XBASE. */
            /*       If XOPT is nonzero, BOBYQB will change it to its usual value later. */
            /*     NF is maintaned as the number of calls of CALFUN so far. */
            /*     KOPT will be such that the least calculated value of F so far is at */
            /*       the point XPT(KOPT,.)+XBASE in the space of the variables. */

            /*     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ, */
            /*     BMAT and ZMAT for the first iteration, and it maintains the values of */
            /*     NF and KOPT. The vector X is also changed by PRELIM. */

            /*     Set some constants. */

            /* Parameter adjustments */
            zmat_dim1 = npt;
            zmat_offset = 1 + zmat_dim1;
            zmat -= zmat_offset;
            xpt_dim1 = npt;
            xpt_offset = 1 + xpt_dim1;
            xpt -= xpt_offset;
            --x;
            --xl;
            --xu;
            --xbase;
            --fval;
            --gopt;
            --hq;
            --pq;
            bmat_dim1 = ndim;
            bmat_offset = 1 + bmat_dim1;
            bmat -= bmat_offset;
            --sl;
            --su;

            /* Function Body */
            half = .5;
            one = 1.;
            two = 2.;
            zero = 0.;
            rhosq = rhobeg * rhobeg;
            recip = one / rhosq;
            np = n + 1;

            /*     Set XBASE to the initial vector of variables, and set the initial */
            /*     elements of XPT, BMAT, HQ, PQ and ZMAT to zero. */

            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                xbase[j] = x[j];
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    /* L10: */
                    xpt[k + j * xpt_dim1] = zero;
                }
                i__2 = ndim;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    /* L20: */
                    bmat[i__ + j * bmat_dim1] = zero;
                }
            }
            i__2 = n * np / 2;
            for (ih = 1; ih <= i__2; ++ih) {
                /* L30: */
                hq[ih] = zero;
            }
            i__2 = npt;
            for (k = 1; k <= i__2; ++k) {
                pq[k] = zero;
                i__1 = npt - np;
                for (j = 1; j <= i__1; ++j) {
                    /* L40: */
                    zmat[k + j * zmat_dim1] = zero;
                }
            }

            /*     Begin the initialization procedure. NF becomes one more than the number */
            /*     of function values so far. The coordinates of the displacement of the */
            /*     next initial interpolation point from XBASE are set in XPT(NF+1,.). */

            nf = 0;
L50:
            nfm = nf;
            nfx = nf - n;
            ++(nf);
            if (nfm <= n << 1) {
                if (nfm >= 1 && nfm <= n) {
                    stepa = rhobeg;
                    if (su[nfm] == zero) {
                        stepa = -stepa;
                    }
                    xpt[nf + nfm * xpt_dim1] = stepa;
                } else if (nfm > n) {
                    stepa = xpt[nf - n + nfx * xpt_dim1];
                    stepb = -(rhobeg);
                    if (sl[nfx] == zero) {
                        /* Computing MIN */
                        d__1 = two * rhobeg, d__2 = su[nfx];
                        stepb = std::min(d__1,d__2);
                    }
                    if (su[nfx] == zero) {
                        /* Computing MAX */
                        d__1 = -two * rhobeg, d__2 = sl[nfx];
                        stepb = std::max(d__1,d__2);
                    }
                    xpt[nf + nfx * xpt_dim1] = stepb;
                }
            } else {
                itemp = (nfm - np) / n;
                jpt = nfm - itemp * n - n;
                ipt = jpt + itemp;
                if (ipt > n) {
                    itemp = jpt;
                    jpt = ipt - n;
                    ipt = itemp;
                }
                xpt[nf + ipt * xpt_dim1] = xpt[ipt + 1 + ipt * xpt_dim1];
                xpt[nf + jpt * xpt_dim1] = xpt[jpt + 1 + jpt * xpt_dim1];
            }

            /*     Calculate the next value of F. The least function value so far and */
            /*     its index are required. */

            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                /* Computing MIN */
                /* Computing MAX */
                d__3 = xl[j], d__4 = xbase[j] + xpt[nf + j * xpt_dim1];
                d__1 = std::max(d__3,d__4), d__2 = xu[j];
                x[j] = std::min(d__1,d__2);
                if (xpt[nf + j * xpt_dim1] == sl[j]) {
                    x[j] = xl[j];
                }
                if (xpt[nf + j * xpt_dim1] == su[j]) {
                    x[j] = xu[j];
                }
                /* L60: */
            }
            f = calfun(mat(&x[1],n));
            fval[nf] = f;
            if (nf == 1) {
                fbeg = f;
                kopt = 1;
            } else if (f < fval[kopt]) {
                kopt = nf;
            }

            /*     Set the nonzero initial elements of BMAT and the quadratic model in the */
            /*     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions */
            /*     of the NF-th and (NF-N)-th interpolation points may be switched, in */
            /*     order that the function value at the first of them contributes to the */
            /*     off-diagonal second derivative terms of the initial quadratic model. */

            if (nf <= (n << 1) + 1) {
                if (nf >= 2 && nf <= n + 1) {
                    gopt[nfm] = (f - fbeg) / stepa;
                    if (npt < nf + n) {
                        bmat[nfm * bmat_dim1 + 1] = -one / stepa;
                        bmat[nf + nfm * bmat_dim1] = one / stepa;
                        bmat[npt + nfm + nfm * bmat_dim1] = -half * rhosq;
                    }
                } else if (nf >= n + 2) {
                    ih = nfx * (nfx + 1) / 2;
                    temp = (f - fbeg) / stepb;
                    diff = stepb - stepa;
                    hq[ih] = two * (temp - gopt[nfx]) / diff;
                    gopt[nfx] = (gopt[nfx] * stepb - temp * stepa) / diff;
                    if (stepa * stepb < zero) {
                        if (f < fval[nf - n]) {
                            fval[nf] = fval[nf - n];
                            fval[nf - n] = f;
                            if (kopt == nf) {
                                kopt = nf - n;
                            }
                            xpt[nf - n + nfx * xpt_dim1] = stepb;
                            xpt[nf + nfx * xpt_dim1] = stepa;
                        }
                    }
                    bmat[nfx * bmat_dim1 + 1] = -(stepa + stepb) / (stepa * stepb);
                    bmat[nf + nfx * bmat_dim1] = -half / xpt[nf - n + nfx *
                        xpt_dim1];
                    bmat[nf - n + nfx * bmat_dim1] = -bmat[nfx * bmat_dim1 + 1] -
                        bmat[nf + nfx * bmat_dim1];
                    zmat[nfx * zmat_dim1 + 1] = std::sqrt(two) / (stepa * stepb);
                    zmat[nf + nfx * zmat_dim1] = std::sqrt(half) / rhosq;
                    zmat[nf - n + nfx * zmat_dim1] = -zmat[nfx * zmat_dim1 + 1] -
                        zmat[nf + nfx * zmat_dim1];
                }

                /*     Set the off-diagonal second derivatives of the Lagrange functions and */
                /*     the initial quadratic model. */

            } else {
                ih = ipt * (ipt - 1) / 2 + jpt;
                zmat[nfx * zmat_dim1 + 1] = recip;
                zmat[nf + nfx * zmat_dim1] = recip;
                zmat[ipt + 1 + nfx * zmat_dim1] = -recip;
                zmat[jpt + 1 + nfx * zmat_dim1] = -recip;
                temp = xpt[nf + ipt * xpt_dim1] * xpt[nf + jpt * xpt_dim1];
                hq[ih] = (fbeg - fval[ipt + 1] - fval[jpt + 1] + f) / temp;
            }
            if (nf < npt && nf < maxfun) {
                goto L50;
            }

        } /* prelim_ */

    // ----------------------------------------------------------------------------------------

        template <typename funct>
        void rescue_ (
            const funct& calfun,
            const integer n,
            const integer npt,
            const doublereal *xl,
            const doublereal *xu,
            const integer maxfun,
            doublereal *xbase,
            doublereal *xpt,
            doublereal *fval,
            doublereal *xopt,
            doublereal *gopt,
            doublereal *hq,
            doublereal *pq,
            doublereal *bmat,
            doublereal *zmat,
            const integer ndim,
            doublereal *sl,
            doublereal *su,
            integer& nf,
            const doublereal delta,
            integer& kopt,
            doublereal *vlag,
            doublereal * ptsaux,
            doublereal *ptsid,
            doublereal *w
        ) const
        {
            /* System generated locals */
            integer xpt_dim1, xpt_offset, bmat_dim1, bmat_offset, zmat_dim1,
            zmat_offset, i__1, i__2, i__3;
            doublereal d__1, d__2, d__3, d__4;


            /* Local variables */
            doublereal f;
            integer i__, j, k, ih, jp, ip, iq, np, iw;
            doublereal xp = 0, xq = 0, den = 0;
            integer ihp = 0;
            doublereal one;
            integer ihq, jpn, kpt;
            doublereal sum = 0, diff = 0, half = 0, beta = 0;
            integer kold;
            doublereal winc;
            integer nrem, knew;
            doublereal temp, bsum;
            integer nptm;
            doublereal zero = 0, hdiag = 0, fbase = 0, sfrac = 0, denom = 0, vquad = 0, sumpq = 0;
            doublereal dsqmin, distsq, vlmxsq;



            /*     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT, */
            /*       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as */
            /*       the corresponding arguments of BOBYQB on the entry to RESCUE. */
            /*     NF is maintained as the number of calls of CALFUN so far, except that */
            /*       NF is set to -1 if the value of MAXFUN prevents further progress. */
            /*     KOPT is maintained so that FVAL(KOPT) is the least calculated function */
            /*       value. Its correct value must be given on entry. It is updated if a */
            /*       new least function value is found, but the corresponding changes to */
            /*       XOPT and GOPT have to be made later by the calling program. */
            /*     DELTA is the current trust region radius. */
            /*     VLAG is a working space vector that will be used for the values of the */
            /*       provisional Lagrange functions at each of the interpolation points. */
            /*       They are part of a product that requires VLAG to be of length NDIM. */
            /*     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and */
            /*       PTSAUX(2,J) specify the two positions of provisional interpolation */
            /*       points when a nonzero step is taken along e_J (the J-th coordinate */
            /*       direction) through XBASE+XOPT, as specified below. Usually these */
            /*       steps have length DELTA, but other lengths are chosen if necessary */
            /*       in order to satisfy the given bounds on the variables. */
            /*     PTSID is also a working space array. It has NPT components that denote */
            /*       provisional new positions of the original interpolation points, in */
            /*       case changes are needed to restore the linear independence of the */
            /*       interpolation conditions. The K-th point is a candidate for change */
            /*       if and only if PTSID(K) is nonzero. In this case let p and q be the */
            /*       integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p */
            /*       and q are both positive, the step from XBASE+XOPT to the new K-th */
            /*       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise */
            /*       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or */
            /*       p=0, respectively. */
            /*     The first NDIM+NPT elements of the array W are used for working space. */
            /*     The final elements of BMAT and ZMAT are set in a well-conditioned way */
            /*       to the values that are appropriate for the new interpolation points. */
            /*     The elements of GOPT, HQ and PQ are also revised to the values that are */
            /*       appropriate to the final quadratic model. */

            /*     Set some constants. */

            /* Parameter adjustments */
            zmat_dim1 = npt;
            zmat_offset = 1 + zmat_dim1;
            zmat -= zmat_offset;
            xpt_dim1 = npt;
            xpt_offset = 1 + xpt_dim1;
            xpt -= xpt_offset;
            --xl;
            --xu;
            --xbase;
            --fval;
            --xopt;
            --gopt;
            --hq;
            --pq;
            bmat_dim1 = ndim;
            bmat_offset = 1 + bmat_dim1;
            bmat -= bmat_offset;
            --sl;
            --su;
            --vlag;
            ptsaux -= 3;
            --ptsid;
            --w;

            /* Function Body */
            half = .5;
            one = 1.;
            zero = 0.;
            np = n + 1;
            sfrac = half / (doublereal) np;
            nptm = npt - np;

            /*     Shift the interpolation points so that XOPT becomes the origin, and set */
            /*     the elements of ZMAT to zero. The value of SUMPQ is required in the */
            /*     updating of HQ below. The squares of the distances from XOPT to the */
            /*     other interpolation points are set at the end of W. Increments of WINC */
            /*     may be added later to these squares to balance the consideration of */
            /*     the choice of point that is going to become current. */

            sumpq = zero;
            winc = zero;
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                distsq = zero;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    xpt[k + j * xpt_dim1] -= xopt[j];
                    /* L10: */
                    /* Computing 2nd power */
                    d__1 = xpt[k + j * xpt_dim1];
                    distsq += d__1 * d__1;
                }
                sumpq += pq[k];
                w[ndim + k] = distsq;
                winc = std::max(winc,distsq);
                i__2 = nptm;
                for (j = 1; j <= i__2; ++j) {
                    /* L20: */
                    zmat[k + j * zmat_dim1] = zero;
                }
            }

            /*     Update HQ so that HQ and PQ define the second derivatives of the model */
            /*     after XBASE has been shifted to the trust region centre. */

            ih = 0;
            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                w[j] = half * sumpq * xopt[j];
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    /* L30: */
                    w[j] += pq[k] * xpt[k + j * xpt_dim1];
                }
                i__1 = j;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    ++ih;
                    /* L40: */
                    hq[ih] = hq[ih] + w[i__] * xopt[j] + w[j] * xopt[i__];
                }
            }

            /*     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and */
            /*     also set the elements of PTSAUX. */

            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                xbase[j] += xopt[j];
                sl[j] -= xopt[j];
                su[j] -= xopt[j];
                xopt[j] = zero;
                /* Computing MIN */
                d__1 = delta, d__2 = su[j];
                ptsaux[(j << 1) + 1] = std::min(d__1,d__2);
                /* Computing MAX */
                d__1 = -(delta), d__2 = sl[j];
                ptsaux[(j << 1) + 2] = std::max(d__1,d__2);
                if (ptsaux[(j << 1) + 1] + ptsaux[(j << 1) + 2] < zero) {
                    temp = ptsaux[(j << 1) + 1];
                    ptsaux[(j << 1) + 1] = ptsaux[(j << 1) + 2];
                    ptsaux[(j << 1) + 2] = temp;
                }
                if ((d__2 = ptsaux[(j << 1) + 2], std::abs(d__2)) < half * (d__1 = ptsaux[(
                            j << 1) + 1], std::abs(d__1))) {
                    ptsaux[(j << 1) + 2] = half * ptsaux[(j << 1) + 1];
                }
                i__2 = ndim;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    /* L50: */
                    bmat[i__ + j * bmat_dim1] = zero;
                }
            }
            fbase = fval[kopt];

            /*     Set the identifiers of the artificial interpolation points that are */
            /*     along a coordinate direction from XOPT, and set the corresponding */
            /*     nonzero elements of BMAT and ZMAT. */

            ptsid[1] = sfrac;
            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                jp = j + 1;
                jpn = jp + n;
                ptsid[jp] = (doublereal) j + sfrac;
                if (jpn <= npt) {
                    ptsid[jpn] = (doublereal) j / (doublereal) np + sfrac;
                    temp = one / (ptsaux[(j << 1) + 1] - ptsaux[(j << 1) + 2]);
                    bmat[jp + j * bmat_dim1] = -temp + one / ptsaux[(j << 1) + 1];
                    bmat[jpn + j * bmat_dim1] = temp + one / ptsaux[(j << 1) + 2];
                    bmat[j * bmat_dim1 + 1] = -bmat[jp + j * bmat_dim1] - bmat[jpn +
                        j * bmat_dim1];
                    zmat[j * zmat_dim1 + 1] = std::sqrt(2.) / (d__1 = ptsaux[(j << 1) + 1]
                                                          * ptsaux[(j << 1) + 2], std::abs(d__1));
                    zmat[jp + j * zmat_dim1] = zmat[j * zmat_dim1 + 1] * ptsaux[(j <<
                                                                                 1) + 2] * temp;
                    zmat[jpn + j * zmat_dim1] = -zmat[j * zmat_dim1 + 1] * ptsaux[(j
                                                                                   << 1) + 1] * temp;
                } else {
                    bmat[j * bmat_dim1 + 1] = -one / ptsaux[(j << 1) + 1];
                    bmat[jp + j * bmat_dim1] = one / ptsaux[(j << 1) + 1];
                    /* Computing 2nd power */
                    d__1 = ptsaux[(j << 1) + 1];
                    bmat[j + npt + j * bmat_dim1] = -half * (d__1 * d__1);
                }
                /* L60: */
            }

            /*     Set any remaining identifiers with their nonzero elements of ZMAT. */

            if (npt >= n + np) {
                i__2 = npt;
                for (k = np << 1; k <= i__2; ++k) {
                    iw = (integer) (((doublereal) (k - np) - half) / (doublereal) (n)
                    );
                    ip = k - np - iw * n;
                    iq = ip + iw;
                    if (iq > n) {
                        iq -= n;
                    }
                    ptsid[k] = (doublereal) ip + (doublereal) iq / (doublereal) np +
                        sfrac;
                    temp = one / (ptsaux[(ip << 1) + 1] * ptsaux[(iq << 1) + 1]);
                    zmat[(k - np) * zmat_dim1 + 1] = temp;
                    zmat[ip + 1 + (k - np) * zmat_dim1] = -temp;
                    zmat[iq + 1 + (k - np) * zmat_dim1] = -temp;
                    /* L70: */
                    zmat[k + (k - np) * zmat_dim1] = temp;
                }
            }
            nrem = npt;
            kold = 1;
            knew = kopt;

            /*     Reorder the provisional points in the way that exchanges PTSID(KOLD) */
            /*     with PTSID(KNEW). */

L80:
            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                temp = bmat[kold + j * bmat_dim1];
                bmat[kold + j * bmat_dim1] = bmat[knew + j * bmat_dim1];
                /* L90: */
                bmat[knew + j * bmat_dim1] = temp;
            }
            i__2 = nptm;
            for (j = 1; j <= i__2; ++j) {
                temp = zmat[kold + j * zmat_dim1];
                zmat[kold + j * zmat_dim1] = zmat[knew + j * zmat_dim1];
                /* L100: */
                zmat[knew + j * zmat_dim1] = temp;
            }
            ptsid[kold] = ptsid[knew];
            ptsid[knew] = zero;
            w[ndim + knew] = zero;
            --nrem;
            if (knew != kopt) {
                temp = vlag[kold];
                vlag[kold] = vlag[knew];
                vlag[knew] = temp;

                /*     Update the BMAT and ZMAT matrices so that the status of the KNEW-th */
                /*     interpolation point can be changed from provisional to original. The */
                /*     branch to label 350 occurs if all the original points are reinstated. */
                /*     The nonnegative values of W(NDIM+K) are required in the search below. */

                update_(n, npt, &bmat[bmat_offset], &zmat[zmat_offset], ndim, &vlag[1],
                        beta, denom, knew, &w[1]);
                if (nrem == 0) {
                    goto L350;
                }
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    /* L110: */
                    w[ndim + k] = (d__1 = w[ndim + k], std::abs(d__1));
                }
            }

            /*     Pick the index KNEW of an original interpolation point that has not */
            /*     yet replaced one of the provisional interpolation points, giving */
            /*     attention to the closeness to XOPT and to previous tries with KNEW. */

L120:
            dsqmin = zero;
            i__2 = npt;
            for (k = 1; k <= i__2; ++k) {
                if (w[ndim + k] > zero) {
                    if (dsqmin == zero || w[ndim + k] < dsqmin) {
                        knew = k;
                        dsqmin = w[ndim + k];
                    }
                }
                /* L130: */
            }
            if (dsqmin == zero) {
                goto L260;
            }

            /*     Form the W-vector of the chosen original interpolation point. */

            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                /* L140: */
                w[npt + j] = xpt[knew + j * xpt_dim1];
            }
            i__2 = npt;
            for (k = 1; k <= i__2; ++k) {
                sum = zero;
                if (k == kopt) {
                } else if (ptsid[k] == zero) {
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j) {
                        /* L150: */
                        sum += w[npt + j] * xpt[k + j * xpt_dim1];
                    }
                } else {
                    ip = (integer) ptsid[k];
                    if (ip > 0) {
                        sum = w[npt + ip] * ptsaux[(ip << 1) + 1];
                    }
                    iq = (integer) ((doublereal) np * ptsid[k] - (doublereal) (ip *
                                                                               np));
                    if (iq > 0) {
                        iw = 1;
                        if (ip == 0) {
                            iw = 2;
                        }
                        sum += w[npt + iq] * ptsaux[iw + (iq << 1)];
                    }
                }
                /* L160: */
                w[k] = half * sum * sum;
            }

            /*     Calculate VLAG and BETA for the required updating of the H matrix if */
            /*     XPT(KNEW,.) is reinstated in the set of interpolation points. */

            i__2 = npt;
            for (k = 1; k <= i__2; ++k) {
                sum = zero;
                i__1 = n;
                for (j = 1; j <= i__1; ++j) {
                    /* L170: */
                    sum += bmat[k + j * bmat_dim1] * w[npt + j];
                }
                /* L180: */
                vlag[k] = sum;
            }
            beta = zero;
            i__2 = nptm;
            for (j = 1; j <= i__2; ++j) {
                sum = zero;
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    /* L190: */
                    sum += zmat[k + j * zmat_dim1] * w[k];
                }
                beta -= sum * sum;
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    /* L200: */
                    vlag[k] += sum * zmat[k + j * zmat_dim1];
                }
            }
            bsum = zero;
            distsq = zero;
            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                sum = zero;
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    /* L210: */
                    sum += bmat[k + j * bmat_dim1] * w[k];
                }
                jp = j + npt;
                bsum += sum * w[jp];
                i__2 = ndim;
                for (ip = npt + 1; ip <= i__2; ++ip) {
                    /* L220: */
                    sum += bmat[ip + j * bmat_dim1] * w[ip];
                }
                bsum += sum * w[jp];
                vlag[jp] = sum;
                /* L230: */
                /* Computing 2nd power */
                d__1 = xpt[knew + j * xpt_dim1];
                distsq += d__1 * d__1;
            }
            beta = half * distsq * distsq + beta - bsum;
            vlag[kopt] += one;

            /*     KOLD is set to the index of the provisional interpolation point that is */
            /*     going to be deleted to make way for the KNEW-th original interpolation */
            /*     point. The choice of KOLD is governed by the avoidance of a small value */
            /*     of the denominator in the updating calculation of UPDATE. */

            denom = zero;
            vlmxsq = zero;
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                if (ptsid[k] != zero) {
                    hdiag = zero;
                    i__2 = nptm;
                    for (j = 1; j <= i__2; ++j) {
                        /* L240: */
                        /* Computing 2nd power */
                        d__1 = zmat[k + j * zmat_dim1];
                        hdiag += d__1 * d__1;
                    }
                    /* Computing 2nd power */
                    d__1 = vlag[k];
                    den = beta * hdiag + d__1 * d__1;
                    if (den > denom) {
                        kold = k;
                        denom = den;
                    }
                }
                /* L250: */
                /* Computing MAX */
                /* Computing 2nd power */
                d__3 = vlag[k];
                d__1 = vlmxsq, d__2 = d__3 * d__3;
                vlmxsq = std::max(d__1,d__2);
            }
            if (denom <= vlmxsq * .01) {
                w[ndim + knew] = -w[ndim + knew] - winc;
                goto L120;
            }
            goto L80;

            /*     When label 260 is reached, all the final positions of the interpolation */
            /*     points have been chosen although any changes have not been included yet */
            /*     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart */
            /*     from the shift of XBASE, the updating of the quadratic model remains to */
            /*     be done. The following cycle through the new interpolation points begins */
            /*     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero, */
            /*     except that a RETURN occurs if MAXFUN prohibits another value of F. */

L260:
            i__1 = npt;
            for (kpt = 1; kpt <= i__1; ++kpt) {
                if (ptsid[kpt] == zero) {
                    goto L340;
                }
                if (nf >= maxfun) {
                    nf = -1;
                    goto L350;
                }
                ih = 0;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    w[j] = xpt[kpt + j * xpt_dim1];
                    xpt[kpt + j * xpt_dim1] = zero;
                    temp = pq[kpt] * w[j];
                    i__3 = j;
                    for (i__ = 1; i__ <= i__3; ++i__) {
                        ++ih;
                        /* L270: */
                        hq[ih] += temp * w[i__];
                    }
                }
                pq[kpt] = zero;
                ip = (integer) ptsid[kpt];
                iq = (integer) ((doublereal) np * ptsid[kpt] - (doublereal) (ip * np))
                    ;
                if (ip > 0) {
                    xp = ptsaux[(ip << 1) + 1];
                    xpt[kpt + ip * xpt_dim1] = xp;
                }
                if (iq > 0) {
                    xq = ptsaux[(iq << 1) + 1];
                    if (ip == 0) {
                        xq = ptsaux[(iq << 1) + 2];
                    }
                    xpt[kpt + iq * xpt_dim1] = xq;
                }

                /*     Set VQUAD to the value of the current model at the new point. */

                vquad = fbase;
                if (ip > 0) {
                    ihp = (ip + ip * ip) / 2;
                    vquad += xp * (gopt[ip] + half * xp * hq[ihp]);
                }
                if (iq > 0) {
                    ihq = (iq + iq * iq) / 2;
                    vquad += xq * (gopt[iq] + half * xq * hq[ihq]);
                    if (ip > 0) {
                        iw = std::max(ihp,ihq) - (i__3 = ip - iq, std::abs(i__3));
                        vquad += xp * xq * hq[iw];
                    }
                }
                i__3 = npt;
                for (k = 1; k <= i__3; ++k) {
                    temp = zero;
                    if (ip > 0) {
                        temp += xp * xpt[k + ip * xpt_dim1];
                    }
                    if (iq > 0) {
                        temp += xq * xpt[k + iq * xpt_dim1];
                    }
                    /* L280: */
                    vquad += half * pq[k] * temp * temp;
                }

                /*     Calculate F at the new interpolation point, and set DIFF to the factor */
                /*     that is going to multiply the KPT-th Lagrange function when the model */
                /*     is updated to provide interpolation to the new function value. */

                i__3 = n;
                for (i__ = 1; i__ <= i__3; ++i__) {
                    /* Computing MIN */
                    /* Computing MAX */
                    d__3 = xl[i__], d__4 = xbase[i__] + xpt[kpt + i__ * xpt_dim1];
                    d__1 = std::max(d__3,d__4), d__2 = xu[i__];
                    w[i__] = std::min(d__1,d__2);
                    if (xpt[kpt + i__ * xpt_dim1] == sl[i__]) {
                        w[i__] = xl[i__];
                    }
                    if (xpt[kpt + i__ * xpt_dim1] == su[i__]) {
                        w[i__] = xu[i__];
                    }
                    /* L290: */
                }
                ++(nf);
                f = calfun(mat(&w[1],n));
                fval[kpt] = f;
                if (f < fval[kopt]) {
                    kopt = kpt;
                }
                diff = f - vquad;

                /*     Update the quadratic model. The RETURN from the subroutine occurs when */
                /*     all the new interpolation points are included in the model. */

                i__3 = n;
                for (i__ = 1; i__ <= i__3; ++i__) {
                    /* L310: */
                    gopt[i__] += diff * bmat[kpt + i__ * bmat_dim1];
                }
                i__3 = npt;
                for (k = 1; k <= i__3; ++k) {
                    sum = zero;
                    i__2 = nptm;
                    for (j = 1; j <= i__2; ++j) {
                        /* L320: */
                        sum += zmat[k + j * zmat_dim1] * zmat[kpt + j * zmat_dim1];
                    }
                    temp = diff * sum;
                    if (ptsid[k] == zero) {
                        pq[k] += temp;
                    } else {
                        ip = (integer) ptsid[k];
                        iq = (integer) ((doublereal) np * ptsid[k] - (doublereal) (ip
                                                                                   * np));
                        ihq = (iq * iq + iq) / 2;
                        if (ip == 0) {
                            /* Computing 2nd power */
                            d__1 = ptsaux[(iq << 1) + 2];
                            hq[ihq] += temp * (d__1 * d__1);
                        } else {
                            ihp = (ip * ip + ip) / 2;
                            /* Computing 2nd power */
                            d__1 = ptsaux[(ip << 1) + 1];
                            hq[ihp] += temp * (d__1 * d__1);
                            if (iq > 0) {
                                /* Computing 2nd power */
                                d__1 = ptsaux[(iq << 1) + 1];
                                hq[ihq] += temp * (d__1 * d__1);
                                iw = std::max(ihp,ihq) - (i__2 = iq - ip, std::abs(i__2));
                                hq[iw] += temp * ptsaux[(ip << 1) + 1] * ptsaux[(iq <<
                                                                                 1) + 1];
                            }
                        }
                    }
                    /* L330: */
                }
                ptsid[kpt] = zero;
L340:
                ;
            }
L350:
            ;
        } /* rescue_ */

    // ----------------------------------------------------------------------------------------

        void trsbox_(
            const integer n,
            const integer npt,
            const doublereal *xpt,
            const doublereal *xopt,
            const doublereal *gopt,
            const doublereal *hq,
            const doublereal *pq,
            const doublereal *sl,
            const doublereal *su,
            const doublereal delta,
            doublereal *xnew,
            doublereal *d__,
            doublereal *gnew,
            doublereal *xbdi,
            doublereal *s,
            doublereal *hs,
            doublereal *hred,
            doublereal *dsq,
            doublereal *crvmin
        ) const
        {
            /* System generated locals */
            integer xpt_dim1, xpt_offset, i__1, i__2;
            doublereal d__1, d__2, d__3, d__4;

            /* Local variables */
            integer i__, j, k, ih;
            doublereal ds;
            integer iu;
            doublereal dhd, dhs, cth, one, shs, sth, ssq, half, beta, sdec, blen;
            integer iact = 0, nact = 0;
            doublereal angt, qred;
            integer isav;
            doublereal temp = 0, zero = 0, xsav = 0, xsum = 0, angbd = 0, dredg = 0, sredg = 0;
            integer iterc;
            doublereal resid = 0, delsq = 0, ggsav = 0, tempa = 0, tempb = 0,
                       redmax = 0, dredsq = 0, redsav = 0, onemin = 0, gredsq = 0, rednew = 0;
            integer itcsav = 0;
            doublereal rdprev = 0, rdnext = 0, stplen = 0, stepsq = 0;
            integer itermax = 0;


            /*     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same */
            /*       meanings as the corresponding arguments of BOBYQB. */
            /*     DELTA is the trust region radius for the present calculation, which */
            /*       seeks a small value of the quadratic model within distance DELTA of */
            /*       XOPT subject to the bounds on the variables. */
            /*     XNEW will be set to a new vector of variables that is approximately */
            /*       the one that minimizes the quadratic model within the trust region */
            /*       subject to the SL and SU constraints on the variables. It satisfies */
            /*       as equations the bounds that become active during the calculation. */
            /*     D is the calculated trial step from XOPT, generated iteratively from an */
            /*       initial value of zero. Thus XNEW is XOPT+D after the final iteration. */
            /*     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated */
            /*       when D is updated. */
            /*     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is */
            /*       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the */
            /*       I-th variable has become fixed at a bound, the bound being SL(I) or */
            /*       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This */
            /*       information is accumulated during the construction of XNEW. */
            /*     The arrays S, HS and HRED are also used for working space. They hold the */
            /*       current search direction, and the changes in the gradient of Q along S */
            /*       and the reduced D, respectively, where the reduced D is the same as D, */
            /*       except that the components of the fixed variables are zero. */
            /*     DSQ will be set to the square of the length of XNEW-XOPT. */
            /*     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise */
            /*       it is set to the least curvature of H that occurs in the conjugate */
            /*       gradient searches that are not restricted by any constraints. The */
            /*       value CRVMIN=-1.0D0 is set, however, if all of these searches are */
            /*       constrained. */

            /*     A version of the truncated conjugate gradient is applied. If a line */
            /*     search is restricted by a constraint, then the procedure is restarted, */
            /*     the values of the variables that are at their bounds being fixed. If */
            /*     the trust region boundary is reached, then further changes may be made */
            /*     to D, each one being in the two dimensional space that is spanned */
            /*     by the current D and the gradient of Q at XOPT+D, staying on the trust */
            /*     region boundary. Termination occurs when the reduction in Q seems to */
            /*     be close to the greatest reduction that can be achieved. */

            /*     Set some constants. */

            /* Parameter adjustments */
            xpt_dim1 = npt;
            xpt_offset = 1 + xpt_dim1;
            xpt -= xpt_offset;
            --xopt;
            --gopt;
            --hq;
            --pq;
            --sl;
            --su;
            --xnew;
            --d__;
            --gnew;
            --xbdi;
            --s;
            --hs;
            --hred;

            /* Function Body */
            half = .5;
            one = 1.;
            onemin = -1.;
            zero = 0.;

            /*     The sign of GOPT(I) gives the sign of the change to the I-th variable */
            /*     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether */
            /*     or not to fix the I-th variable at one of its bounds initially, with */
            /*     NACT being set to the number of fixed variables. D and GNEW are also */
            /*     set for the first iteration. DELSQ is the upper bound on the sum of */
            /*     squares of the free variables. QRED is the reduction in Q so far. */

            iterc = 0;
            nact = 0;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                xbdi[i__] = zero;
                if (xopt[i__] <= sl[i__]) {
                    if (gopt[i__] >= zero) {
                        xbdi[i__] = onemin;
                    }
                } else if (xopt[i__] >= su[i__]) {
                    if (gopt[i__] <= zero) {
                        xbdi[i__] = one;
                    }
                }
                if (xbdi[i__] != zero) {
                    ++nact;
                }
                d__[i__] = zero;
                /* L10: */
                gnew[i__] = gopt[i__];
            }
            delsq = delta * delta;
            qred = zero;
            *crvmin = onemin;

            /*     Set the next search direction of the conjugate gradient method. It is */
            /*     the steepest descent direction initially and when the iterations are */
            /*     restarted because a variable has just been fixed by a bound, and of */
            /*     course the components of the fixed variables are zero. ITERMAX is an */
            /*     upper bound on the indices of the conjugate gradient iterations. */

L20:
            beta = zero;
L30:
            stepsq = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (xbdi[i__] != zero) {
                    s[i__] = zero;
                } else if (beta == zero) {
                    s[i__] = -gnew[i__];
                } else {
                    s[i__] = beta * s[i__] - gnew[i__];
                }
                /* L40: */
                /* Computing 2nd power */
                d__1 = s[i__];
                stepsq += d__1 * d__1;
            }
            if (stepsq == zero) {
                goto L190;
            }
            if (beta == zero) {
                gredsq = stepsq;
                itermax = iterc + n - nact;
            }
            if (gredsq * delsq <= qred * 1e-4 * qred) {
                goto L190;
            }

            /*     Multiply the search direction by the second derivative matrix of Q and */
            /*     calculate some scalars for the choice of steplength. Then set BLEN to */
            /*     the length of the the step to the trust region boundary and STPLEN to */
            /*     the steplength, ignoring the simple bounds. */

            goto L210;
L50:
            resid = delsq;
            ds = zero;
            shs = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (xbdi[i__] == zero) {
                    /* Computing 2nd power */
                    d__1 = d__[i__];
                    resid -= d__1 * d__1;
                    ds += s[i__] * d__[i__];
                    shs += s[i__] * hs[i__];
                }
                /* L60: */
            }
            if (resid <= zero) {
                goto L90;
            }
            temp = std::sqrt(stepsq * resid + ds * ds);
            if (ds < zero) {
                blen = (temp - ds) / stepsq;
            } else {
                blen = resid / (temp + ds);
            }
            stplen = blen;
            if (shs > zero) {
                /* Computing MIN */
                d__1 = blen, d__2 = gredsq / shs;
                stplen = std::min(d__1,d__2);
            }

            /*     Reduce STPLEN if necessary in order to preserve the simple bounds, */
            /*     letting IACT be the index of the new constrained variable. */

            iact = 0;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (s[i__] != zero) {
                    xsum = xopt[i__] + d__[i__];
                    if (s[i__] > zero) {
                        temp = (su[i__] - xsum) / s[i__];
                    } else {
                        temp = (sl[i__] - xsum) / s[i__];
                    }
                    if (temp < stplen) {
                        stplen = temp;
                        iact = i__;
                    }
                }
                /* L70: */
            }

            /*     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q. */

            sdec = zero;
            if (stplen > zero) {
                ++iterc;
                temp = shs / stepsq;
                if (iact == 0 && temp > zero) {
                    *crvmin = std::min(*crvmin,temp);
                    if (*crvmin == onemin) {
                        *crvmin = temp;
                    }
                }
                ggsav = gredsq;
                gredsq = zero;
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    gnew[i__] += stplen * hs[i__];
                    if (xbdi[i__] == zero) {
                        /* Computing 2nd power */
                        d__1 = gnew[i__];
                        gredsq += d__1 * d__1;
                    }
                    /* L80: */
                    d__[i__] += stplen * s[i__];
                }
                /* Computing MAX */
                d__1 = stplen * (ggsav - half * stplen * shs);
                sdec = std::max(d__1,zero);
                qred += sdec;
            }

            /*     Restart the conjugate gradient method if it has hit a new bound. */

            if (iact > 0) {
                ++nact;
                xbdi[iact] = one;
                if (s[iact] < zero) {
                    xbdi[iact] = onemin;
                }
                /* Computing 2nd power */
                d__1 = d__[iact];
                delsq -= d__1 * d__1;
                if (delsq <= zero) {
                    goto L90;
                }
                goto L20;
            }

            /*     If STPLEN is less than BLEN, then either apply another conjugate */
            /*     gradient iteration or RETURN. */

            if (stplen < blen) {
                if (iterc == itermax) {
                    goto L190;
                }
                if (sdec <= qred * .01) {
                    goto L190;
                }
                beta = gredsq / ggsav;
                goto L30;
            }
L90:
            *crvmin = zero;

            /*     Prepare for the alternative iteration by calculating some scalars */
            /*     and by multiplying the reduced D by the second derivative matrix of */
            /*     Q, where S holds the reduced D in the call of GGMULT. */

L100:
            if (nact >= n - 1) {
                goto L190;
            }
            dredsq = zero;
            dredg = zero;
            gredsq = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (xbdi[i__] == zero) {
                    /* Computing 2nd power */
                    d__1 = d__[i__];
                    dredsq += d__1 * d__1;
                    dredg += d__[i__] * gnew[i__];
                    /* Computing 2nd power */
                    d__1 = gnew[i__];
                    gredsq += d__1 * d__1;
                    s[i__] = d__[i__];
                } else {
                    s[i__] = zero;
                }
                /* L110: */
            }
            itcsav = iterc;
            goto L210;

            /*     Let the search direction S be a linear combination of the reduced D */
            /*     and the reduced G that is orthogonal to the reduced D. */

L120:
            ++iterc;
            temp = gredsq * dredsq - dredg * dredg;
            if (temp <= qred * 1e-4 * qred) {
                goto L190;
            }
            temp = std::sqrt(temp);
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (xbdi[i__] == zero) {
                    s[i__] = (dredg * d__[i__] - dredsq * gnew[i__]) / temp;
                } else {
                    s[i__] = zero;
                }
                /* L130: */
            }
            sredg = -temp;

            /*     By considering the simple bounds on the variables, calculate an upper */
            /*     bound on the tangent of half the angle of the alternative iteration, */
            /*     namely ANGBD, except that, if already a free variable has reached a */
            /*     bound, there is a branch back to label 100 after fixing that variable. */

            angbd = one;
            iact = 0;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (xbdi[i__] == zero) {
                    tempa = xopt[i__] + d__[i__] - sl[i__];
                    tempb = su[i__] - xopt[i__] - d__[i__];
                    if (tempa <= zero) {
                        ++nact;
                        xbdi[i__] = onemin;
                        goto L100;
                    } else if (tempb <= zero) {
                        ++nact;
                        xbdi[i__] = one;
                        goto L100;
                    }
                    /* Computing 2nd power */
                    d__1 = d__[i__];
                    /* Computing 2nd power */
                    d__2 = s[i__];
                    ssq = d__1 * d__1 + d__2 * d__2;
                    /* Computing 2nd power */
                    d__1 = xopt[i__] - sl[i__];
                    temp = ssq - d__1 * d__1;
                    if (temp > zero) {
                        temp = std::sqrt(temp) - s[i__];
                        if (angbd * temp > tempa) {
                            angbd = tempa / temp;
                            iact = i__;
                            xsav = onemin;
                        }
                    }
                    /* Computing 2nd power */
                    d__1 = su[i__] - xopt[i__];
                    temp = ssq - d__1 * d__1;
                    if (temp > zero) {
                        temp = std::sqrt(temp) + s[i__];
                        if (angbd * temp > tempb) {
                            angbd = tempb / temp;
                            iact = i__;
                            xsav = one;
                        }
                    }
                }
                /* L140: */
            }

            /*     Calculate HHD and some curvatures for the alternative iteration. */

            goto L210;
L150:
            shs = zero;
            dhs = zero;
            dhd = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (xbdi[i__] == zero) {
                    shs += s[i__] * hs[i__];
                    dhs += d__[i__] * hs[i__];
                    dhd += d__[i__] * hred[i__];
                }
                /* L160: */
            }

            /*     Seek the greatest reduction in Q for a range of equally spaced values */
            /*     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of */
            /*     the alternative iteration. */

            redmax = zero;
            isav = 0;
            redsav = zero;
            iu = (integer) (angbd * 17. + 3.1);
            i__1 = iu;
            for (i__ = 1; i__ <= i__1; ++i__) {
                angt = angbd * (doublereal) i__ / (doublereal) iu;
                sth = (angt + angt) / (one + angt * angt);
                temp = shs + angt * (angt * dhd - dhs - dhs);
                rednew = sth * (angt * dredg - sredg - half * sth * temp);
                if (rednew > redmax) {
                    redmax = rednew;
                    isav = i__;
                    rdprev = redsav;
                } else if (i__ == isav + 1) {
                    rdnext = rednew;
                }
                /* L170: */
                redsav = rednew;
            }

            /*     Return if the reduction is zero. Otherwise, set the sine and cosine */
            /*     of the angle of the alternative iteration, and calculate SDEC. */

            if (isav == 0) {
                goto L190;
            }
            if (isav < iu) {
                temp = (rdnext - rdprev) / (redmax + redmax - rdprev - rdnext);
                angt = angbd * ((doublereal) isav + half * temp) / (doublereal) iu;
            }
            cth = (one - angt * angt) / (one + angt * angt);
            sth = (angt + angt) / (one + angt * angt);
            temp = shs + angt * (angt * dhd - dhs - dhs);
            sdec = sth * (angt * dredg - sredg - half * sth * temp);
            if (sdec <= zero) {
                goto L190;
            }

            /*     Update GNEW, D and HRED. If the angle of the alternative iteration */
            /*     is restricted by a bound on a free variable, that variable is fixed */
            /*     at the bound. */

            dredg = zero;
            gredsq = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                gnew[i__] = gnew[i__] + (cth - one) * hred[i__] + sth * hs[i__];
                if (xbdi[i__] == zero) {
                    d__[i__] = cth * d__[i__] + sth * s[i__];
                    dredg += d__[i__] * gnew[i__];
                    /* Computing 2nd power */
                    d__1 = gnew[i__];
                    gredsq += d__1 * d__1;
                }
                /* L180: */
                hred[i__] = cth * hred[i__] + sth * hs[i__];
            }
            qred += sdec;
            if (iact > 0 && isav == iu) {
                ++nact;
                xbdi[iact] = xsav;
                goto L100;
            }

            /*     If SDEC is sufficiently small, then RETURN after setting XNEW to */
            /*     XOPT+D, giving careful attention to the bounds. */

            if (sdec > qred * .01) {
                goto L120;
            }
L190:
            *dsq = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                /* Computing MAX */
                /* Computing MIN */
                d__3 = xopt[i__] + d__[i__], d__4 = su[i__];
                d__1 = std::min(d__3,d__4), d__2 = sl[i__];
                xnew[i__] = std::max(d__1,d__2);
                if (xbdi[i__] == onemin) {
                    xnew[i__] = sl[i__];
                }
                if (xbdi[i__] == one) {
                    xnew[i__] = su[i__];
                }
                d__[i__] = xnew[i__] - xopt[i__];
                /* L200: */
                /* Computing 2nd power */
                d__1 = d__[i__];
                *dsq += d__1 * d__1;
            }
            return;
            /*     The following instructions multiply the current S-vector by the second */
            /*     derivative matrix of the quadratic model, putting the product in HS. */
            /*     They are reached from three different parts of the software above and */
            /*     they can be regarded as an external subroutine. */

L210:
            ih = 0;
            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                hs[j] = zero;
                i__2 = j;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    ++ih;
                    if (i__ < j) {
                        hs[j] += hq[ih] * s[i__];
                    }
                    /* L220: */
                    hs[i__] += hq[ih] * s[j];
                }
            }
            i__2 = npt;
            for (k = 1; k <= i__2; ++k) {
                if (pq[k] != zero) {
                    temp = zero;
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j) {
                        /* L230: */
                        temp += xpt[k + j * xpt_dim1] * s[j];
                    }
                    temp *= pq[k];
                    i__1 = n;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        /* L240: */
                        hs[i__] += temp * xpt[k + i__ * xpt_dim1];
                    }
                }
                /* L250: */
            }
            if (*crvmin != zero) {
                goto L50;
            }
            if (iterc > itcsav) {
                goto L150;
            }
            i__2 = n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                /* L260: */
                hred[i__] = hs[i__];
            }
            goto L120;
        } /* trsbox_ */

    // ----------------------------------------------------------------------------------------

        void update_(
            const integer n,
            const integer npt,
            doublereal *bmat,
            doublereal *zmat,
            const integer ndim,
            doublereal *vlag,
            const doublereal beta,
            const doublereal denom,
            const integer knew,
            doublereal *w
        ) const
        {
            /* System generated locals */
            integer bmat_dim1, bmat_offset, zmat_dim1, zmat_offset, i__1, i__2;
            doublereal d__1, d__2, d__3;

            /* Local variables */
            integer i__, j, k, jp;
            doublereal one, tau, temp;
            integer nptm;
            doublereal zero, alpha, tempa, tempb, ztest;


            /*     The arrays BMAT and ZMAT are updated, as required by the new position */
            /*     of the interpolation point that has the index KNEW. The vector VLAG has */
            /*     N+NPT components, set on entry to the first NPT and last N components */
            /*     of the product Hw in equation (4.11) of the Powell (2006) paper on */
            /*     NEWUOA. Further, BETA is set on entry to the value of the parameter */
            /*     with that name, and DENOM is set to the denominator of the updating */
            /*     formula. Elements of ZMAT may be treated as zero if their moduli are */
            /*     at most ZTEST. The first NDIM elements of W are used for working space. */

            /*     Set some constants. */

            /* Parameter adjustments */
            zmat_dim1 = npt;
            zmat_offset = 1 + zmat_dim1;
            zmat -= zmat_offset;
            bmat_dim1 = ndim;
            bmat_offset = 1 + bmat_dim1;
            bmat -= bmat_offset;
            --vlag;
            --w;

            /* Function Body */
            one = 1.;
            zero = 0.;
            nptm = npt - n - 1;
            ztest = zero;
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                i__2 = nptm;
                for (j = 1; j <= i__2; ++j) {
                    /* L10: */
                    /* Computing MAX */
                    d__2 = ztest, d__3 = (d__1 = zmat[k + j * zmat_dim1], std::abs(d__1));
                    ztest = std::max(d__2,d__3);
                }
            }
            ztest *= 1e-20;

            /*     Apply the rotations that put zeros in the KNEW-th row of ZMAT. */

            i__2 = nptm;
            for (j = 2; j <= i__2; ++j) {
                if ((d__1 = zmat[knew + j * zmat_dim1], std::abs(d__1)) > ztest) {
                    /* Computing 2nd power */
                    d__1 = zmat[knew + zmat_dim1];
                    /* Computing 2nd power */
                    d__2 = zmat[knew + j * zmat_dim1];
                    temp = std::sqrt(d__1 * d__1 + d__2 * d__2);
                    tempa = zmat[knew + zmat_dim1] / temp;
                    tempb = zmat[knew + j * zmat_dim1] / temp;
                    i__1 = npt;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        temp = tempa * zmat[i__ + zmat_dim1] + tempb * zmat[i__ + j *
                            zmat_dim1];
                        zmat[i__ + j * zmat_dim1] = tempa * zmat[i__ + j * zmat_dim1]
                            - tempb * zmat[i__ + zmat_dim1];
                        /* L20: */
                        zmat[i__ + zmat_dim1] = temp;
                    }
                }
                zmat[knew + j * zmat_dim1] = zero;
                /* L30: */
            }

            /*     Put the first NPT components of the KNEW-th column of HLAG into W, */
            /*     and calculate the parameters of the updating formula. */

            i__2 = npt;
            for (i__ = 1; i__ <= i__2; ++i__) {
                w[i__] = zmat[knew + zmat_dim1] * zmat[i__ + zmat_dim1];
                /* L40: */
            }
            alpha = w[knew];
            tau = vlag[knew];
            vlag[knew] -= one;

            /*     Complete the updating of ZMAT. */

            temp = std::sqrt(denom);
            tempb = zmat[knew + zmat_dim1] / temp;
            tempa = tau / temp;
            i__2 = npt;
            for (i__ = 1; i__ <= i__2; ++i__) {
                /* L50: */
                zmat[i__ + zmat_dim1] = tempa * zmat[i__ + zmat_dim1] - tempb * vlag[
                    i__];
            }

            /*     Finally, update the matrix BMAT. */

            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                jp = npt + j;
                w[jp] = bmat[knew + j * bmat_dim1];
                tempa = (alpha * vlag[jp] - tau * w[jp]) / denom;
                tempb = (-(beta) * w[jp] - tau * vlag[jp]) / denom;
                i__1 = jp;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    bmat[i__ + j * bmat_dim1] = bmat[i__ + j * bmat_dim1] + tempa *
                        vlag[i__] + tempb * w[i__];
                    if (i__ > npt) {
                        bmat[jp + (i__ - npt) * bmat_dim1] = bmat[i__ + j *
                            bmat_dim1];
                    }
                    /* L60: */
                }
            }
        } /* update_ */
    };

// ----------------------------------------------------------------------------------------

    template <
        typename funct,
        typename T,
        typename U
        >
    double find_min_bobyqaConvgFail (
        const funct& f,
        T& x,
        long npt,
        const U& x_lower,
        const U& x_upper,
        const double rho_begin,
        const double rho_end,
        const long max_f_evals,  int& ConvgFail
    )
    {
        // The starting point (i.e. x) must be a column vector.
        COMPILE_TIME_ASSERT(T::NC <= 1);

        // check the requirements.  Also split the assert up so that the error message isn't huge.
        DLIB_CASSERT(is_col_vector(x) && is_col_vector(x_lower) && is_col_vector(x_upper) &&
                    x.size() == x_lower.size() && x_lower.size() == x_upper.size() &&
                    x.size() > 1 && max_f_evals > 1,
            "\tvoid find_min_bobyqa()"
            << "\n\t Invalid arguments have been given to this function"
            << "\n\t is_col_vector(x):       " << is_col_vector(x)
            << "\n\t is_col_vector(x_lower): " << is_col_vector(x_lower)
            << "\n\t is_col_vector(x_upper): " << is_col_vector(x_upper)
            << "\n\t x.size():               " << x.size()
            << "\n\t x_lower.size():         " << x_lower.size()
            << "\n\t x_upper.size():         " << x_upper.size()
            << "\n\t max_f_evals:            " << max_f_evals
        );


        DLIB_CASSERT(x.size() + 2 <= npt && npt <= (x.size()+1)*(x.size()+2)/2 &&
                    0 < rho_end && rho_end < rho_begin &&
                    min(x_upper - x_lower) > 2*rho_begin &&
                    min(x - x_lower) >= 0 && min(x_upper - x) >= 0,
            "\tvoid find_min_bobyqa()"
            << "\n\t Invalid arguments have been given to this function"
            << "\n\t ntp in valid range: " << (x.size() + 2 <= npt && npt <= (x.size()+1)*(x.size()+2)/2)
            << "\n\t npt:                " << npt
            << "\n\t rho_begin:          " << rho_begin
            << "\n\t rho_end:            " << rho_end
            << "\n\t min(x_upper - x_lower) > 2*rho_begin:           " << (min(x_upper - x_lower) > 2*rho_begin)
            << "\n\t min(x - x_lower) >= 0 && min(x_upper - x) >= 0: " << (min(x - x_lower) >= 0 && min(x_upper - x) >= 0)
        );
       /*

        if(x.size() + 2 <= npt && npt <= (x.size()+1)*(x.size()+2)/2 &&
                    0 < rho_end && rho_end < rho_begin &&
                    min(x_upper - x_lower) > 2*rho_begin &&
                    min(x - x_lower) >= 0 && min(x_upper - x) >= 0)
        {
          std::cout<<"Invalid arguments given for find_min_bobyqaConvgFail!!!\n";
          std::cout<<"0 < rho_end && rho_end < rho_begin &&\
                    min(x_upper - x_lower) > 2*rho_begin &&\
                    min(x - x_lower) >= 0 && min(x_upper - x) >= 0 fails\n";
        ConvgFail = 1; return 0;
        }
        */


        bobyqa_implementationConvgFail impl;
        return impl.find_min(f, x, npt, x_lower, x_upper, rho_begin, rho_end, max_f_evals,ConvgFail);
    }



  template <
        typename search_strategy_type,
        typename stop_strategy_type,
        typename funct,
        typename funct_der,
        typename T,
        typename EXP1,
        typename EXP2
        >
    double find_min_box_constrainedConvFail (
        search_strategy_type search_strategy,
        stop_strategy_type stop_strategy,
        const funct& f,
        const funct_der& der,
        T& x,
        const matrix_exp<EXP1>& x_lower,
        const matrix_exp<EXP2>& x_upper, int& ConvgFail, const int& MaxIteration
    )
    {
        /*
            The implementation of this function is more or less based on the discussion in
            the paper Projected Newton-type Methods in Machine Learning by Mark Schmidt, et al.
        */

        // make sure the requires clause is not violated
        COMPILE_TIME_ASSERT(is_matrix<T>::value);
        // The starting point (i.e. x) must be a column vector.
        COMPILE_TIME_ASSERT(T::NC <= 1);

        DLIB_CASSERT (
            is_col_vector(x) && is_col_vector(x_lower) && is_col_vector(x_upper) &&
            x.size() == x_lower.size() && x.size() == x_upper.size(),
            "\tdouble find_min_box_constrained()"
            << "\n\t The inputs to this function must be equal length column vectors."
            << "\n\t is_col_vector(x):       " << is_col_vector(x)
            << "\n\t is_col_vector(x_upper): " << is_col_vector(x_upper)
            << "\n\t is_col_vector(x_upper): " << is_col_vector(x_upper)
            << "\n\t x.size():               " << x.size()
            << "\n\t x_lower.size():         " << x_lower.size()
            << "\n\t x_upper.size():         " << x_upper.size()
        );
        DLIB_ASSERT (
            min(x_upper-x_lower) > 0,
            "\tdouble find_min_box_constrained()"
            << "\n\t You have to supply proper box constraints to this function."
            << "\n\r min(x_upper-x_lower): " << min(x_upper-x_lower)
        );


        T g, s;
        double f_value = f(x);
        g = der(x);

        if (!is_finite(f_value))
          // throw error("The objective function generated non-finite outputs");
		{  ConvgFail=1;   return 0; }
        if (!is_finite(g))
         //   throw error("The objective function generated non-finite outputs");
      {  ConvgFail=1;   return 0; }

        // gap_eps determines how close we have to get to a bound constraint before we
        // start basically dropping it from the optimization and consider it to be an
        // active constraint.
        const double gap_eps = 1e-8;

        double last_alpha = 1;  int IterDone=0;
        while(stop_strategy.should_continue_search(x, f_value, g))
        {
if(IterDone>MaxIteration)  {  ConvgFail=1;   return 0; }

            s = search_strategy.get_next_direction(x, f_value, zero_bounded_variables(gap_eps, g, x, g, x_lower, x_upper));
            s = gap_step_assign_bounded_variables(gap_eps, s, x, g, x_lower, x_upper);

            double alpha = backtracking_line_search(
                        make_line_search_function(clamp_function(f,x_lower,x_upper), x, s, f_value),
                        f_value,
                        dot(g,s), // compute gradient for the line search
                        last_alpha,
                        search_strategy.get_wolfe_rho(),
                        search_strategy.get_max_line_search_iterations());

            // Do a trust region style thing for alpha.  The idea is that if we take a
            // small step then we are likely to take another small step.  So we reuse the
            // alpha from the last iteration unless the line search didn't shrink alpha at
            // all, in that case, we start with a bigger alpha next time.
            if (alpha == last_alpha)
                last_alpha = std::min(last_alpha*10,1.0);
            else
                last_alpha = alpha;

            // Take the search step indicated by the above line search
            x = clamp(x + alpha*s, x_lower, x_upper);
            g = der(x);

  if (!is_finite(f_value))
            // throw error("The objective function generated non-finite outputs");
			 { ConvgFail=1; return 0; }
 if (!is_finite(g))
             // throw error("The objective function generated non-finite outputs");
  { ConvgFail=1; return 0; }

	      IterDone=IterDone+1;
        }

        return f_value;
    }


    template <
        typename search_strategy_type,
        typename stop_strategy_type,
        typename funct,
        typename funct_der,
        typename T
        >
    double find_minConvFail(
        search_strategy_type search_strategy,
        stop_strategy_type stop_strategy,
        const funct& f,
        const funct_der& der,
        T& x,
        double min_f, int& ConvFail, const int& MaxIter
    )
    {
        COMPILE_TIME_ASSERT(is_matrix<T>::value);
        // The starting point (i.e. x) must be a column vector.
        COMPILE_TIME_ASSERT(T::NC <= 1);

        DLIB_CASSERT (
            is_col_vector(x),
            "\tdouble find_min()"
            << "\n\tYou have to supply column vectors to this function"
            << "\n\tx.nc():    " << x.nc()
        );

        T g, s;

        double f_value = f(x);
        g = der(x);

        if (!is_finite(f_value))
        // throw error("The objective function generated non-finite outputs");
{  ConvFail=1; return 0;  }
        if (!is_finite(g))
         //  throw error("The objective function generated non-finite outputs");
{  ConvFail=1; return 0;  }

      int IterCounts=0;  int MaxIterPlus1=MaxIter+1;
        while(stop_strategy.should_continue_search(x, f_value, g) && f_value > min_f )
        {
            IterCounts = IterCounts+1;
       if( IterCounts > MaxIterPlus1)  { ConvFail=1; return 0;  }
            s = search_strategy.get_next_direction(x, f_value, g);

            double alpha = line_search(
                        make_line_search_function(f,x,s, f_value),
                        f_value,
                        make_line_search_function(der,x,s, g),
                        dot(g,s), // compute initial gradient for the line search
                        search_strategy.get_wolfe_rho(), search_strategy.get_wolfe_sigma(), min_f,
                        search_strategy.get_max_line_search_iterations());

            // Take the search step indicated by the above line search
            x += alpha*s;

            if (!is_finite(f_value))
             // throw error("The objective function generated non-finite outputs");
{  ConvFail=1; return 0;  }
            if (!is_finite(g))
            //  throw error("The objective function generated non-finite outputs");
{  ConvFail=1; return 0;  }
        }

        return f_value;
    }




   template <
        typename stop_strategy_type,
        typename funct_model
        >
    double find_min_trust_regionConvFail (
        stop_strategy_type stop_strategy,
        const funct_model& model,
        typename funct_model::column_vector& x,\
         int& ConvFail, const int& MaxIteration,\
        double radius = 1
    )
    {
        /*
            This is an implementation of algorithm 4.1(Trust Region)
            from the book Numerical Optimization by Nocedal and Wright.
        */

        // make sure requires clause is not broken
        DLIB_ASSERT(is_col_vector(x) && radius > 0,
            "\t double find_min_trust_region()"
            << "\n\t invalid arguments were given to this function"
            << "\n\t is_col_vector(x): " << is_col_vector(x)
            << "\n\t radius:           " << radius
            );

        const double initial_radius = radius;

        typedef typename funct_model::column_vector T;
        typedef typename T::type type;

        typename funct_model::general_matrix h;
        typename funct_model::column_vector g, p, d;
        type f_value = model(x);

        model.get_derivative_and_hessian(x,g,h);

 // DLIB_ASSERT(is_finite(x), "The objective function generated non-finite outputs");

 if( !is_finite(x) )  {ConvFail=1; return 0; }

 // DLIB_ASSERT(is_finite(g), "The objective function generated non-finite outputs");

 if( !is_finite(g) )  {ConvFail=1; return 0; }

//  DLIB_ASSERT(is_finite(h), "The objective function generated non-finite outputs");

 if( !is_finite(h) )  {ConvFail=1; return 0; }

        // Sometimes the loop below won't modify x because the trust region step failed.
        // This bool tells us when we are in that case.
        bool stale_x = false;
   int IterDone=0;
        while(stale_x || stop_strategy.should_continue_search(x, f_value, g))
        {

     if(IterDone>MaxIteration) {ConvFail=1; return 0; }
            const unsigned long iter = solve_trust_region_subproblem(h,
                                                                     g,
                                                                     radius,
                                                                     p,
                                                                     0.1,
                                                                     20);


            const type new_f_value = model(x+p);
            const type predicted_improvement = -0.5*trans(p)*h*p - trans(g)*p;
            const type measured_improvement = (f_value - new_f_value);

            // If the sub-problem can't find a way to improve then stop.  This only happens when p is essentially 0.
            if (std::abs(predicted_improvement) <= std::abs(measured_improvement)*std::numeric_limits<type>::epsilon())
                break;

            // predicted_improvement shouldn't be negative but it might be if something went
            // wrong in the trust region solver.  So put abs() here to guard against that.  This
            // way the sign of rho is determined only by the sign of measured_improvement.
            const type rho = measured_improvement/std::abs(predicted_improvement);



            if (rho < 0.25)
            {
                radius *= 0.25;

                // something has gone horribly wrong if the radius has shrunk to zero.  So just
                // give up if that happens.
                if (radius <= initial_radius*std::numeric_limits<double>::epsilon())
                    break;
            }
            else
            {
                // if rho > 0.75 and we are being checked by the radius
                if (rho > 0.75 && iter > 1)
                {
                    radius = std::min<type>(1000,  2*radius);
                }
            }

            if (rho > 0)
            {
                x = x + p;
                f_value = new_f_value;
                model.get_derivative_and_hessian(x,g,h);
                stale_x = false;
            }
            else
            {
                stale_x = true;
            }


//  DLIB_ASSERT(is_finite(x), "The objective function generated non-finite outputs");

 if( !is_finite(x) )  {ConvFail=1; return 0; }


//  DLIB_ASSERT(is_finite(g), "The objective function generated non-finite outputs");

if( !is_finite(g) )  {ConvFail=1; return 0; }


// DLIB_ASSERT(is_finite(h), "The objective function generated non-finite outputs");

if( !is_finite(h) )  {ConvFail=1; return 0; }

    IterDone=IterDone+1;

        }

        return f_value;
    }


/*  Compute the machine epsilon.
   The epsilon function returns the machine epsilon: the smallest number that, when summed with 1, produces a value greater than one.
  */
  // alternatively use std::numeric_limits::epsilon
template <typename T>
T epsilonGet()
{
    T eps, del, neweps;
    del    = (T) 0.5;
    eps    = (T) 0.0;
    neweps = (T) 1.0;

    while ( del > 0 )
    {
        if ( 1 + neweps > 1 )    /* Then the value might be too large */
        {
            eps = neweps;    /* ...save the current value... */
            neweps -= del;    /* ...and decrement a bit */
        }
        else          /* Then the value is too small */
        {
            neweps += del;    /* ...so increment it */
        }
        del *= 0.5;      /* Reduce the adjustment by half */
    }
    return eps;
}

  double* p_step_size;

// const double h2 = std::sqrt(epsilonGet<double>());
// const double h = std::sqrt(h2);

// numerical gradient
dbcolvec gradcdif_old(double (*FUNCTOR) (const dbcolvec& ), const dbcolvec& th)
{   double h = *p_step_size;  double f_f, f_b;
    int dim = th.nr();
    dbcolvec grad(dim);
    dbcolvec th_F(dim);
    dbcolvec th_B(dim);
    for(int i = 0; i< dim; ++i)
    {   th_F = th; th_B = th;
        th_F(i) = th(i) + h;  f_f = FUNCTOR(th_F);
        th_B(i) = th(i) - h;  f_b = FUNCTOR(th_B);
        grad(i) = (f_f - f_b) / (2*h);

         }
    return grad;
}


dbcolvec gradcdif(double (*FUNCTOR) (const dbcolvec& ),  const dbcolvec& th)
{   double h;
    h = *p_step_size; // 1e-7 is the DEFAULT of derivative() in dlib
    h=0.001;
double f_f, f_b;
    int dim = th.nr();

    dbcolvec grad(dim);
    dbcolvec th_F(dim);
    dbcolvec th_B(dim);

    for(int i = 0; i< dim; ++i)
    {   th_F = th; th_B = th;
        th_F(i) = th(i) + h;  f_f = FUNCTOR(th_F);
        th_B(i) = th(i) - h;  f_b = FUNCTOR(th_B);
        grad(i) = (f_f - f_b) / (2*h);
  // if(dim==22) { DebugOut<<"\ni = "<<i<<"\n"<<"f(x+ei) = "<<f_f<<"\n"<<"f(x-ei) = "<<f_b<<"\n"<<"grad(i) = "<<grad(i)<<"\n"; }
    }
    return grad;
}


// five points formula at https://en.wikipedia.org/wiki/Five-point_stencil
dbcolvec gradcdif_5p(double (*FUNCTOR) (const dbcolvec& ),  const dbcolvec& th)
{
      double h;
    h = *p_step_size; // 1e-7 is the DEFAULT of derivative() in dlib
    // h=0.0001;
double f_f, f_b, f_f2h, f_b2h;
    int dim = th.nr();

    dbcolvec grad(dim), th_temp(dim);

    for(int i = 0; i< dim; ++i)
    {   th_temp =  th;
        th_temp(i) = th(i) + h;
        f_f = FUNCTOR(th_temp);
        th_temp =  th;
        th_temp(i) = th(i) - h;
        f_b = FUNCTOR(th_temp);
        th_temp =  th;
        th_temp(i) = th(i) + 2*h;
        f_f2h = FUNCTOR(th_temp);
        th_temp =  th;
        th_temp(i) = th(i) - 2*h;
        f_b2h = FUNCTOR(th_temp);
        grad(i) = (8*f_f + f_b2h - 8*f_b -f_f2h) / (12*h);
    if(dim==22) { DebugOut<<"\ni = "<<i<<"\n"<<"f(x+ei) = "<<f_f<<"\n"<<"f(x-ei) = "<<f_b<<"\n"<<"f(x+2ei) = "<<f_f2h<<"\n"<<"f(x-2ei) = "<<f_b2h<<"\n"; }
    }
    return grad;
}

// numerical hessian
dbmat HESScdif(double (*fun) (const dbcolvec& ),  const dbcolvec& th)
{
  // double h = 0.00001;
    double h = *p_step_size;
    double h2 = h*h;
    int dim = th.nr();
    dbmat hess(dim,dim);
    dbcolvec ei(dim);
    dbcolvec ej(dim);
    double fval = fun(th);
    for ( int i = 0; i < dim; ++i)
    {
        ei=0;
        ei(i)=h;
        for ( int j = 0; j < dim; ++j)
        {
            ej = 0;
            ej(j) = h;

            if (i == j)
            {
                hess(i,i) = ( -fun(th + 2.0 * ei) + 16.0 *
                              fun(th + ei) - 30.0 * fval + 16.0 *
                              fun(th - ei) -
                              fun(th - 2.0 * ei)) / (12.0 * h2);
            }
            else
            {
               if(i<j)
                { hess(i,j) = ( fun(th + ei + ej) - fun(th + ei - ej)
                              - fun(th - ei + ej) + fun(th - ei - ej))
                            / (4.0 * h2);
                }
              else  { hess(i,j) = hess(j,i);   }

            }
        }
    }


    return hess;
}


dbmat numerial_hess_using_analytical_grad(const dbcolvec& ci,double (*FUNCTOR) (const dbcolvec& ))
 {
   dbmat temp(4,4);
   return temp;
 }




/*
 Eta is dim*(dim+1)/2 by 1; L is dim*dim
 correspondence between Eta and L (L is such that Sig=L*t(L) ):
  Eta(a_s,0) corresponds to  L(s,s)   then a_s - a_{s-1}= p+1-s; a0=0;
 a_s=dim*s+s-s(s+1)/2;  for example, dim=3; a1=3
 //input: Eta is dim*(dim+1)/2 by 1;
//output: L is dim*dim;
 */
void EtaToL( const dbcolvec& Eta, const int& dim, dbmat& L)
{
    L=0;
    for(int s=0; s<dim; ++s)
    {
        int a_s=dim*s+s-s*(s+1)/2;
        L(s,s)=exp(Eta(a_s));
        if(s<dim-1)
            set_subm(L,range(s+1,dim-1), range(s,s)) = subm(Eta,range(a_s+1,a_s+dim-s-1), range(0,0));
        // L(s+1,s,dim-1,s)=Eta(a_s+1,0,a_s+dim-s-1,0);
        // note a_{s+1} - a_s= dim-s
    }
}

//L to Eta
void LToEta(const dbmat& L, const int& dim, dbcolvec& Eta )
{
    for(int s=0; s<dim; ++s)
    {
        int a_s=dim*s+s-s*(s+1)/2;
        Eta(a_s)=log( L(s,s));
        if(s<dim-1)
            set_subm(Eta,range(a_s+1,a_s+dim-s-1), range(0,0))=subm(L,range(s+1,dim-1), range(s,s));
    }
}

// initialize the matrix for quadrature points
void GetQHPts(std::ifstream& QHfile, int& m, dbmat& QHpts)
  {
    int i,j;
    for(i=0;i<m;++i)
    {
     for(j=0;j<2;++j)
    QHfile>>QHpts(i,j);
    }
  }

void GetMinMax(const dbcolvec& Input, double& Min, double& Max)
{
   int len = Input.nr();  Min = Input(0); Max = Input(0);
   for(int k=0; k<len; ++k)
   {
       if( Input(k) < Min) Min = Input(k);
       else if(Input(k) > Max) Max = Input(k);
   }
}


void GetMinMax_with_index(const dbcolvec& Input, double& Min, double& Max, int& min_ind, int&max_ind)
{
   int len = Input.nr();
    Min = Input(0);
    Max = Input(0);
    min_ind = 0;
    max_ind = 0;
   for(int k=0; k<len; ++k)
   {
       if( Input(k) < Min) { Min = Input(k);  min_ind = k;  }
       else if(Input(k) > Max) { Max = Input(k);   max_ind = k;  }
   }
}



// input: Dat; row i for obs i (1*dim);
// output: sample mean vector mu; sample Covariance matrix Sig
void SamMeanCov(const dbmat& Dat, const int& SamSz, const int& dim,\
    dbcolvec& mu, dbmat& Sig)
{
    for(int i=0; i<dim; ++i) mu(i)=mean( subm(Dat,range(0,SamSz-1),range(i,i)) );
    dbmat All1(SamSz,1);
    All1=1;
    dbmat SamMeanExt(dim,SamSz);
    SamMeanExt=All1*trans(mu);
    Sig=trans(Dat-SamMeanExt)*(Dat-SamMeanExt) /SamSz;
}





 // global pointers
 dbmat* pDat_;  dbmat* p_GH_points;
 dbcolvec* p_ni_vec;
 int* p_nsub;  const dbcolvec* p_para_vec;
 int* p_start_j_sub_i; int* p_endj_sub_i;
 dbmat* p_cov_ranf_inv; double* p_det_cov_ranf_inv;
 dlib::rand* p_rd; dbcolvec* p_std_mvn_sam; dbcolvec* p_ranf_mvn_sam;
 dbcolvec* p_sam_mean; dbmat* p_sam_cov;
 int* p_num_mc; double* p_neg_loglike_U95; double* p_neg_loglike_L95;


 // pDat_ columns:  0-id, 1-u_ij, 2-v_ij, 3-t_ij, 4-HT_i, 5-LT_i, first 1/3 Control group, second 1/3 LT, last 1/3 HT
   // para_vec: 0-3 fixed effects for log-odds-omit; 4-7 fixed effects for amount-omit; 8:log_sd_eij; 9-11 eta determining cov(bi0,bi1)
 void jm_gen_simu_data(dlib::rand&rd, const int&nsub, const dbcolvec&ni_vec, const dbcolvec&para_vec_true, const dbmat&cov_ranef_L, dbmat&dat)
{
   int i, j, s, startj_subi, endj_subi, nsub3 = 3*nsub;
   double tij, linear_pred_censor,   HTi, LTi, positive_prob, unif_rand;
   dbcolvec ranef_realization(2);

   startj_subi = 0;
   for(i=0; i<nsub3; ++i)
  {
    if(i>0) {  startj_subi = startj_subi + int(ni_vec(i-1)); }
    endj_subi = startj_subi + int(ni_vec(i));


    if(i<nsub)
   {
     HTi = 0; LTi = 0;
   }

    if(i>(nsub-1) & i<(2*nsub) )
   {
     HTi = 0; LTi = 1;
   }

    if( i>(2*nsub-1) & i<(3*nsub)  )
   {
     HTi = 1; LTi = 0;
   }

     // generate realizations of random effects
     for(s=0;s<2;++s)
     { ranef_realization(s) = rd.get_random_gaussian(); }

      ranef_realization = cov_ranef_L*ranef_realization;

   for(j=startj_subi; j<endj_subi; ++j)
   {
       dat(j,0) = i+1;
       tij = double(j-startj_subi)/ni_vec(i);
       dat(j,3) = tij;
   linear_pred_censor = ranef_realization(0) +  para_vec_true(0) + para_vec_true(1)*tij + para_vec_true(2)*HTi + para_vec_true(3)*LTi;
   positive_prob = 1/(  1+std::exp(-1.0*linear_pred_censor) );
    unif_rand = rd.get_random_double();

     if(unif_rand < positive_prob)  { dat(j,1) = 1;
     dat(j,2) = ranef_realization(1) +  para_vec_true(4) + para_vec_true(5)*tij + para_vec_true(6)*HTi + para_vec_true(7)*LTi + std::exp(para_vec_true(8))*rd.get_random_gaussian();

     }
      else  {
            dat(j,1) = 0;
      dat(j,2) = 0;
      }

     dat(j,4) = HTi;
     dat(j,5) = LTi;
    }
  }

}

 // pDat_ columns:  0-id, 1-u_ij, 2-v_ij, 3-t_ij, 4-HT_i, 5-LT_i, first 1/3 Control group, second 1/3 LT, last 1/3 HT
   // para_vec: 0-3 fixed effects for log-odds-omit; 4-7 fixed effects for amount-omit; 8:log_sd_eij; 9-11 eta determining cov(bi0,bi1)
 double MC_neg_loglike(const dbcolvec& para_vec)
 {
     double loglike, integral_i, loglike_i, num_mc = double(*p_num_mc), condi_loglikei, temp;
     double  residual, m2_i, se, linear_pred_censored, tij, HTi, LTi;
     int start_j_sub_i, ni, k, j, endj ,s;
     dbcolvec cov_ranef_Lvec(3);  dbmat cov_ranef_L(2,2), cov_ranef_inv(2,2);
     cov_ranef_Lvec = subm(para_vec,range(9,11),range(0,0));
     int dim_rf = 2; EtaToL(cov_ranef_Lvec,dim_rf,cov_ranef_L);
//  std::cout<<"cov_ranef_L\n"<<cov_ranef_L;

      loglike = 0;
      start_j_sub_i = 0;
      *p_neg_loglike_U95 = 0;
      *p_neg_loglike_L95 = 0;
       *p_sam_mean = 0; *p_sam_cov = 0;
    for(int i = 0; i < 1; ++i) // *p_nsub
    {
         if(i > 0)
        { start_j_sub_i = start_j_sub_i \
         + int( (*p_ni_vec)(i-1) );  }
         ni = int( (*p_ni_vec)(i) );
         endj = start_j_sub_i + ni;
 start_j_sub_i = 0; endj = 1;
         integral_i = 0; m2_i = 0;
      for(k=0; k<(*p_num_mc);++k)
      {
       for(s=0; s<dim_rf; ++s)
       { (*p_std_mvn_sam)(s) = p_rd->get_random_gaussian(); }
       *p_ranf_mvn_sam = cov_ranef_L*(*p_std_mvn_sam);
        condi_loglikei = 0;

        for(j = start_j_sub_i; j < endj; ++j)//
        {
             tij = (*pDat_)(j,3);
             HTi = (*pDat_)(j,4);
             LTi = (*pDat_)(j,5);

    // if( k==0) std::cout<<"tij="<<tij<<"; j="<<j<<"; (*pDat_)(j,1)="<<(*pDat_)(j,1)<<"\n";
            linear_pred_censored = (*p_ranf_mvn_sam)(0) +  para_vec(0) + para_vec(1)*tij + para_vec(2)*HTi + para_vec(3)*LTi;
            if (std::abs((*pDat_)(j,1)) <0.1)
             {
                condi_loglikei = condi_loglikei \
               - std::log(1 + std::exp(linear_pred_censored));
             }
             else
             {
                 condi_loglikei = condi_loglikei \
                 + linear_pred_censored - std::log(1 + std::exp(linear_pred_censored));
              residual = (*pDat_)(j,2) -  para_vec(4) - para_vec(5)*tij - para_vec(6)*HTi - para_vec(7)*LTi -(*p_ranf_mvn_sam)(1);
              condi_loglikei = condi_loglikei -para_vec(8) - 0.5*residual*residual/(std::exp(2*para_vec(8) ));
          }
        }

        integral_i = integral_i + std::exp(condi_loglikei);
        m2_i = m2_i + std::exp(condi_loglikei*2);
       *p_sam_mean = *p_sam_mean + *p_ranf_mvn_sam;
      *p_sam_cov = *p_sam_cov + (*p_ranf_mvn_sam)*trans(*p_ranf_mvn_sam);
   //if(k<5){ std::cout<<"k="<<k<<";  exp(condi_loglikei)="<<std::exp(condi_loglikei)<<"; exp(condi_loglikei*condi_loglikei)="<<std::exp(condi_loglikei*condi_loglikei)<<"\n";  }
      }
     *p_sam_mean = *p_sam_mean/num_mc;
     *p_sam_cov = *p_sam_cov/num_mc - (*p_sam_mean)*trans(*p_sam_mean);
       integral_i = integral_i/num_mc;
       m2_i = m2_i/num_mc;
  //std::cout<<"i="<<i<<"; integral_i="<<integral_i<<"\n";
  //std::cout<<"i="<<i<<"; m2_i="<<m2_i<<"\n";
       se = std::sqrt(   (m2_i-integral_i*integral_i)/num_mc);
       loglike_i = std::log(integral_i);
 // std::cout<<"i="<<i<<"; MC_loglike_i="<<loglike_i<<"\n";
 // std::cout<<"i="<<i<<"; condi_loglikei="<<(condi_loglikei)<<"\n";
 //std::cout<<"i="<<i<<"; var="<<(m2_i-integral_i*integral_i)<<"\n";
       loglike = loglike + loglike_i;
  // std::cout<<"MC, i="<<i<<"; loglike="<<loglike_i<<"\n";
       *p_neg_loglike_U95 = *p_neg_loglike_U95 + std::log(integral_i+1.96*se);
       *p_neg_loglike_L95 = *p_neg_loglike_L95 + std::log(integral_i-1.96*se);
    }
     temp = *p_neg_loglike_U95;
    *p_neg_loglike_U95 = (-1.0)*(*p_neg_loglike_L95);
    *p_neg_loglike_L95 = (-1.0)*(temp);
 return (-1.0)*loglike;
 }

  int GH_neg_loglike_direct_d2_counter;
  double GH_neg_loglike_direct_d2(const dbcolvec& para_vec)
  {
     double loglike, integral_i, loglike_i, condi_loglikei;
     double  residual, linear_pred_censored, tij, HTi, LTi;
     int start_j_sub_i, ni, j, endj;

     int dim_rf = 2;
     double half_dim_rf = double(dim_rf)/2;

     dbcolvec cov_ranef_Lvec(3);
     dbmat cov_ranef_L(2,2), cov_ranef_inv(2,2);
     cov_ranef_Lvec = subm(para_vec,range(9,11),range(0,0));
     EtaToL(cov_ranef_Lvec,dim_rf,cov_ranef_L);
//  std::cout<<"inside GH neg loglike, cov_ranef_L\n"<<cov_ranef_L<<"\n";
//   std::cout<<"(*p_GH_points)\n"<<(*p_GH_points);

     int num_GH_pts = (*p_GH_points).nr(), s0, s1;
     dbcolvec z_star(dim_rf), z_shift(dim_rf);

     loglike = 0;
     start_j_sub_i = 0;
    for(int i = 0; i < 1; ++i)  // i < *p_nsub;
    {
         if(i > 0)
        { start_j_sub_i = start_j_sub_i \
         + int( (*p_ni_vec)(i-1) );  }
         ni = int( (*p_ni_vec)(i) );
         endj = start_j_sub_i + ni;
       //endj = 5;
    start_j_sub_i = 0;  endj = 1;
    integral_i = 0;
    for(s0 = 0; s0 < num_GH_pts; ++s0)
    {   z_star(0) = (*p_GH_points)(s0,0);

      for(s1 = 0; s1 < num_GH_pts; ++s1)
    {   z_star(1) = (*p_GH_points)(s1,0);

        z_shift = std::sqrt(2)*cov_ranef_L*z_star;

       condi_loglikei = 0;
        for(j = start_j_sub_i; j < endj; ++j) //
        {
             tij = (*pDat_)(j,3);
             HTi = (*pDat_)(j,4);
             LTi = (*pDat_)(j,5);

  // std::cout<<"GH, tij="<<tij<<"; j="<<j<<"; (*pDat_)(j,1)="<<(*pDat_)(j,1)<<"\n";
            linear_pred_censored = para_vec(0) + para_vec(1)*tij + para_vec(2)*HTi + para_vec(3)*LTi + z_shift(0) ;
            if (std::abs((*pDat_)(j,1)) <0.1)
             {
                condi_loglikei = condi_loglikei \
                 - std::log(1 + std::exp(linear_pred_censored));
             }
             else
             {
              condi_loglikei = condi_loglikei \
               + linear_pred_censored - std::log(1 + std::exp(linear_pred_censored));
              residual = (*pDat_)(j,2) -  para_vec(4) - para_vec(5)*tij - para_vec(6)*HTi - para_vec(7)*LTi - z_shift(1);
             condi_loglikei = condi_loglikei  -para_vec(8)- 0.5*residual*residual/(std::exp(2*para_vec(8) ));
          }
        }

    integral_i = integral_i + std::exp(condi_loglikei)*((*p_GH_points)(s0,1))\
  *((*p_GH_points)(s1,1));
    }
    }


    integral_i = integral_i*std::pow(dlib::pi,-1.0*half_dim_rf);
  // std::cout<<"GH, i="<<i<<"; integral_i="<<integral_i<<"\n";
       loglike_i =  std::log(integral_i);
 // std::cout<<"GH, i="<<i<<"; loglike="<<loglike_i<<"\n";
       loglike = loglike + loglike_i;
    }
    GH_neg_loglike_direct_d2_counter++;
   //   if( (GH_neg_loglike_direct_d2_counter%100==0) ) {
   // std::cout<<"#like_eva="<<GH_neg_loglike_direct_d2_counter<<"; loglike="<<(-1.0)*loglike<<"; parameter0-5 ="<<trans(subm(para_vec,range(0,5),range(0,0)) ); }
 // DebugOut<<"loglike = "<<loglike<<"\n";
 // if(!is_finite(loglike)) { DebugOut<<std::setprecision(12)<<"para_vec with infinity likelihood\n"<<trans(para_vec);  }
 // std::cout<<"loglike = "<<loglike<<"\n";
 return (-1.0)*loglike;
  }


 dbcolvec Grad_GH_neg_loglike_direct_d2(const dbcolvec& para_vec)
  { return gradcdif(GH_neg_loglike_direct_d2,para_vec);}




  // for prob of vomiting: para_vec 0-3 fixed effects para_vec are random effects
  double neg_loglike_glmm_d1(const dbcolvec&para_vec)
 {
     double loglike, integral_i, loglike_i, condi_loglikei;
     double   linear_pred_censored, tij, HTi, LTi;
     int start_j_sub_i, ni, j, endj;

     int dim_rf = 1;
     double half_dim_rf = double(dim_rf)/2;

     dbcolvec cov_ranef_Lvec(1);
     dbmat cov_ranef_L(1,1), cov_ranef_inv(1,1);
     cov_ranef_Lvec = subm(para_vec,range(4,4),range(0,0));
     EtaToL(cov_ranef_Lvec,dim_rf,cov_ranef_L);
//  std::cout<<"inside GH neg loglike, cov_ranef_L\n"<<cov_ranef_L<<"\n";
//   std::cout<<"(*p_GH_points)\n"<<(*p_GH_points);

     int num_GH_pts = (*p_GH_points).nr(), s0;
     dbcolvec z_star(dim_rf), z_shift(dim_rf);

     loglike = 0;
     start_j_sub_i = 0;

      for(int i = 0; i < *p_nsub; ++i)  // i < *p_nsub;
    {
         if(i > 0)
        { start_j_sub_i = start_j_sub_i \
         + int( (*p_ni_vec)(i-1) );  }
         ni = int( (*p_ni_vec)(i) );
         endj = start_j_sub_i + ni;
       //endj = 5;
  // start_j_sub_i = 0;  endj = 7;
    integral_i = 0;
    for(s0 = 0; s0 < num_GH_pts; ++s0)
    {   z_star(0) = (*p_GH_points)(s0,0);

        z_shift = std::sqrt(2)*cov_ranef_L*z_star;

       condi_loglikei = 0;
        for(j = start_j_sub_i; j < endj; ++j) //
        {
             tij = (*pDat_)(j,3);
             HTi = (*pDat_)(j,4);
             LTi = (*pDat_)(j,5);

  // std::cout<<"GH, tij="<<tij<<"; j="<<j<<"; (*pDat_)(j,1)="<<(*pDat_)(j,1)<<"\n";
            linear_pred_censored = para_vec(0) + para_vec(1)*tij + para_vec(2)*HTi + para_vec(3)*LTi + z_shift(0) ;
            if (std::abs((*pDat_)(j,1)) <0.1)
             {
                condi_loglikei = condi_loglikei \
                 - std::log(1 + std::exp(linear_pred_censored));
             }
             else
             {
              condi_loglikei = condi_loglikei \
            + linear_pred_censored - std::log(1 + std::exp(linear_pred_censored));

          }
        }

    integral_i = integral_i + std::exp(condi_loglikei)*((*p_GH_points)(s0,1));
    }



    integral_i = integral_i*std::pow(dlib::pi,-1.0*half_dim_rf);
  // std::cout<<"GH, i="<<i<<"; integral_i="<<integral_i<<"\n";
       loglike_i =  std::log(integral_i);
 // std::cout<<"GH, i="<<i<<"; loglike="<<loglike_i<<"\n";
       loglike = loglike + loglike_i;
    }
 return (-1.0)*loglike;
 }


  // for prob of vomiting:
  //  para_vec 0-1 fixed effects intercept and slope for tij
  // para_vec_2: coef_HTi-coef_LTi, para_vec_3 Coef LTi
  double neg_loglike_glmm_d1_ord(const dbcolvec&para_vec)
 {
     double loglike, integral_i, loglike_i, condi_loglikei;
     double  linear_pred_censored, tij, HTi, LTi;
     int start_j_sub_i, ni, j, endj;

     int dim_rf = 1;
     double half_dim_rf = double(dim_rf)/2;
     double coef_HTi = para_vec(2) + para_vec(3), coef_LTi = para_vec(3);

     dbcolvec cov_ranef_Lvec(1);
     dbmat cov_ranef_L(1,1), cov_ranef_inv(1,1);
     cov_ranef_Lvec = subm(para_vec,range(4,4),range(0,0));
     EtaToL(cov_ranef_Lvec,dim_rf,cov_ranef_L);
//  std::cout<<"inside GH neg loglike, cov_ranef_L\n"<<cov_ranef_L<<"\n";
//   std::cout<<"(*p_GH_points)\n"<<(*p_GH_points);

     int num_GH_pts = (*p_GH_points).nr(), s0;
     dbcolvec z_star(dim_rf), z_shift(dim_rf);

     loglike = 0;
     start_j_sub_i = 0;

      for(int i = 0; i < *p_nsub; ++i)  // i < *p_nsub;
    {
         if(i > 0)
        { start_j_sub_i = start_j_sub_i \
         + int( (*p_ni_vec)(i-1) );  }
         ni = int( (*p_ni_vec)(i) );
         endj = start_j_sub_i + ni;
       //endj = 5;
   // start_j_sub_i = 0;  endj = 1;
    integral_i = 0;
    for(s0 = 0; s0 < num_GH_pts; ++s0)
    {   z_star(0) = (*p_GH_points)(s0,0);

        z_shift = std::sqrt(2)*cov_ranef_L*z_star;

       condi_loglikei = 0;
        for(j = start_j_sub_i; j < endj; ++j) //
        {
             tij = (*pDat_)(j,3);
             HTi = (*pDat_)(j,4);
             LTi = (*pDat_)(j,5);


            linear_pred_censored = para_vec(0) + para_vec(1)*tij + coef_HTi*HTi + coef_LTi*LTi + z_shift(0) ;

     //  if(s0<1) {  std::cout<<"GH GLMM, s0 = "<<s0<<";  coef_HTi = "<< coef_HTi<<";  coef_LTi = "<< coef_LTi<<"\n"<<"para_vec(0) = "<<para_vec(0)<<"; para_vec(1) = " <<para_vec(1)<<"; GH GLMM, s0 = "<<s0<<";  tij = "<< tij<<";  HTi = "<< HTi<<";  LTi = "<< LTi<<"\n"<<"s0 = "<<s0<<";  z_shift(0) = "<< z_shift(0)<<"; GH GLMM, linear_pred_censored = "<<linear_pred_censored<<"\n";  }

            if (std::abs((*pDat_)(j,1)) <0.1)
             {
                condi_loglikei = condi_loglikei \
                 - std::log(1 + std::exp(linear_pred_censored));
             }
             else
             {
              condi_loglikei = condi_loglikei \
            + linear_pred_censored - std::log(1 + std::exp(linear_pred_censored));

          }
        }
    // std::cout<<"s0 = "<<s0<<"; condi_loglikei = "<<condi_loglikei<<"\n";
    integral_i = integral_i + std::exp(condi_loglikei)*((*p_GH_points)(s0,1));
    }



    integral_i = integral_i*std::pow(dlib::pi,-1.0*half_dim_rf);
  // std::cout<<"GH, i="<<i<<"; integral_i="<<integral_i<<"\n";
       loglike_i =  std::log(integral_i);
 // std::cout<<"GH, i="<<i<<"; loglike="<<loglike_i<<"\n";
       loglike = loglike + loglike_i;
    }
 return (-1.0)*loglike;
 }




  // for linear amount volume:   //  para_vec 0-1 fixed effects intercept and slope for tij
  // para_vec_2: coef_HTi-coef_LTi, para_vec_3 Coef LTi
  // para_vec_4 are random effects sd, para_5 are sd_eij;
  // analytical LME intergral:
  // key  \int exp(-Ax^2-Bx-C) dx = exp(B^2/(4*A) - C)*\sqrt{1/(2A)} sqrt{2pi}
   double neg_loglike_lme_ana_ord(const dbcolvec& para_vec)
    {
     double loglike, loglike_i, sig_b = std::exp(para_vec(4)), mu_b, sig_e = std::exp( para_vec(5) );
     int start_j_sub_i, ni, j, endj, ni_star;
     double coef_A, coef_B, coef_C;
     double vij_minus_xij_beta, tij, HTi, LTi;
     double coef_HTi = para_vec(2) + para_vec(3), coef_LTi = para_vec(3);

     mu_b = 0;

     int dim_rf = 1;

     dbcolvec cov_ranef_Lvec(1);
     dbmat cov_ranef_L(1,1), cov_ranef_inv(1,1);
     cov_ranef_Lvec = subm(para_vec,range(4,4),range(0,0));
     EtaToL(cov_ranef_Lvec,dim_rf,cov_ranef_L);

     loglike = 0;
     start_j_sub_i = 0;
       for(int i = 0;i < *p_nsub; ++i)  // i < *p_nsub;
    {
         if(i > 0)
        { start_j_sub_i = start_j_sub_i \
         + int( (*p_ni_vec)(i-1) );  }

         ni = int( (*p_ni_vec)(i) );
    // std::cout<<"ni = "<<ni<<"\n";
         endj = start_j_sub_i + ni;
  // start_j_sub_i = 0; endj = 1;

       coef_A = 0; coef_B = 0; coef_C = 0; ni_star = 0;

        for(j = start_j_sub_i; j < endj; ++j) //
        {


       if (std::abs( (*pDat_)(j,1) ) < 0.1)
             {
                continue;
             }
          else
          {
             ni_star = ni_star + 1;
             tij = (*pDat_)(j,3);
             HTi = (*pDat_)(j,4);
             LTi = (*pDat_)(j,5);
   //std::cout<<"Ana, j="<<j<<"; (*pDat_)(j,2)="<<(*pDat_)(j,2)<<"\n";
          vij_minus_xij_beta =  (*pDat_)(j,2) - para_vec(0) - para_vec(1)*tij -coef_HTi*HTi - coef_LTi*LTi;
       coef_A = coef_A + 1/(2*sig_e*sig_e);
       coef_B = coef_B -  vij_minus_xij_beta/(sig_e*sig_e);
       coef_C = coef_C + vij_minus_xij_beta*vij_minus_xij_beta/(2*sig_e*sig_e);
          }
        }
        coef_A = coef_A + 1/(2*sig_b*sig_b);
        coef_B = coef_B - mu_b/(sig_b*sig_b);
        coef_C = coef_C + mu_b*mu_b/(2*sig_b*sig_b);

       loglike_i = -0.5*double(ni_star)*std::log(2*dlib::pi) - ni_star*para_vec(5) + coef_B*coef_B/(4*coef_A)-coef_C + 0.5*std::log(2*dlib::pi) - 0.5*std::log(2*coef_A)-0.5*std::log(2*dlib::pi)-para_vec(4);
     //std::cout<<"Ana, i="<<i<<"; loglike_i="<<loglike_i<<"\n";
       loglike = loglike + loglike_i;
    }
 return (-1.0)*loglike;
 }


   double neg_loglike_lme_ana_glmm_GH_ord(const dbcolvec& para_vec)
    {
     double loglike, weighted_sum;
     double  tij, HTi, LTi;
     int start_j_sub_i, ni, j, endj;
     double coef_HTi_o = para_vec(2) + para_vec(3), coef_LTi_o = para_vec(3);
     double coef_HTi_v = para_vec(6) + para_vec(7), coef_LTi_v = para_vec(7);
  // std::cout<<"coef_HTi_o = "<<coef_HTi_o<<"; "<<"coef_LTi_o = "<<coef_LTi_o<<"; "<<"coef_HTi_v = "<<coef_HTi_v<<"; "<<"coef_LTi_v = "<<coef_LTi_v<<"\n";

     int dim_rf = 2, ni_star;

     dbcolvec cov_ranef_Lvec(3);
     dbmat cov_ranef_L(2,2), cov_ranef(2,2);
     cov_ranef_Lvec = subm(para_vec,range(9,11),range(0,0));
     EtaToL(cov_ranef_Lvec,dim_rf,cov_ranef_L);
  //std::cout<<"inside GH neg neg_loglike_lme_ana_glmm_GH_ord, cov_ranef_L\n"<<cov_ranef_L<<"\n";

     cov_ranef = cov_ranef_L*trans(cov_ranef_L);
 // std::cout<<"inside GH neg neg_loglike_lme_ana_glmm_GH_ord, cov_ranef\n"<<cov_ranef<<"\n";
     double sig_0_sqare = cov_ranef(0,0), sig_1_sqare = cov_ranef(1,1), sig_01 = cov_ranef(0,1);
     double sig_0 = std::sqrt(cov_ranef(0,0));
     double sig_e = std::exp(para_vec(8));
  //  std::cout<<"sig_0 = "<<sig_0<<"\n";
 //   std::cout<<"sig_e = "<<sig_e<<"\n";

     int num_GH_pts = (*p_GH_points).nr(), s0;
  // std::cout<<"(*p_GH_points)\n"<<(*p_GH_points);
  // std::cout<<"num_GH_pts = "<<num_GH_pts<<"\n";
     dbcolvec z_star(1), z_shift(1);


     double coef_A, coef_B, coef_C, u1_c, sig_1c = std::sqrt(sig_1_sqare-sig_01*sig_01/sig_0_sqare), log_sig_1c = std::log(sig_1c);
     double log_prod_f_uij_given_bi0, log_ana_integrated_bi0, linear_pred_vomit, vij_minus_xij_beta;
   //  std::cout<<"sig_1c = "<<sig_1c<<"\n";


     loglike = 0;
     start_j_sub_i = 0;
       for(int i = 0; i < *p_nsub; ++i)  // i < *p_nsub;
    {
         if(i > 0)
        { start_j_sub_i = start_j_sub_i \
         + int( (*p_ni_vec)(i-1) );  }

         ni = int( (*p_ni_vec)(i) );
    // std::cout<<"ni = "<<ni<<"\n";
         endj = start_j_sub_i + ni;
 // start_j_sub_i = 0; endj = 1;

     weighted_sum = 0;
     for(s0 = 0; s0 < num_GH_pts; ++s0)
    {
        z_star(0) = (*p_GH_points)(s0,0);
        z_shift(0) = std::sqrt(2)*sig_0*z_star(0);

       u1_c = sig_01*z_shift(0)/sig_0_sqare;
  //  std::cout<<"u1_c = "<<u1_c<<"\n";

       log_prod_f_uij_given_bi0 = 0;
       coef_A = 0; coef_B = 0; coef_C = 0; ni_star = 0;
        for(j = start_j_sub_i; j < endj; ++j) //
        {
             tij = (*pDat_)(j,3);
             HTi = (*pDat_)(j,4);
             LTi = (*pDat_)(j,5);
         linear_pred_vomit = z_shift(0) + para_vec(0) + para_vec(1)*tij + coef_HTi_o*HTi + coef_LTi_o*LTi;

     // if(s0<1) {    std::cout<<"GH joint, s0 = "<<s0<<";  coef_HTi_o = "<< coef_HTi_o<<";  coef_LTi_o = "<< coef_LTi_o<<"\n"<<"para_vec(0) = "<<para_vec(0)<<"; para_vec(1) = " <<para_vec(1)<<"; GH joint, s0 = "<<s0<<";  tij = "<< tij<<";  HTi = "<< HTi<<";  LTi = "<< LTi<<"\n"; std::cout<<"s0 = "<<s0<<";  z_shift(0) = "<< z_shift(0)<<"; GH joint, linear_pred_vomit = "<<linear_pred_vomit<<"\n"; }

       if (std::abs( (*pDat_)(j,1) ) < 0.1)
             {
       log_prod_f_uij_given_bi0 = log_prod_f_uij_given_bi0 -1.0*std::log(1+std::exp(linear_pred_vomit));
    // std::cout<<" (*pDat_)(j,1) = "<<(*pDat_)(j,1)<<"\n";
             }
          else
          {
       log_prod_f_uij_given_bi0 = log_prod_f_uij_given_bi0 +linear_pred_vomit-1.0*std::log(1+std::exp(linear_pred_vomit));
       ni_star = ni_star + 1;

   // std::cout<<"ni_star, j="<<j<<"; ni_star="<<ni_star<<"\n";
      vij_minus_xij_beta =  (*pDat_)(j,2) - para_vec(4) - para_vec(5)*tij -coef_HTi_v*HTi - coef_LTi_v*LTi;
       coef_A = coef_A + 1/(2*sig_e*sig_e);
       coef_B = coef_B -  vij_minus_xij_beta/(sig_e*sig_e);
       coef_C = coef_C + vij_minus_xij_beta*vij_minus_xij_beta/(2*sig_e*sig_e);
          }
        }

        coef_A = coef_A + 1/(2*sig_1c*sig_1c);
        coef_B = coef_B - u1_c/(sig_1c*sig_1c);
        coef_C = coef_C + u1_c*u1_c/(2*sig_1c*sig_1c);
   //  std::cout<<"coef_A = "<<coef_A<<"; coef_B = "<<coef_B<<"; coef_C = "<<coef_C<<"\n";

        log_ana_integrated_bi0 = -0.5*double(ni_star)*std::log(2*dlib::pi) - double(ni_star)*std::log(sig_e) + coef_B*coef_B/(4*coef_A)-coef_C + 0.5*std::log(2*dlib::pi) - 0.5*std::log(2*coef_A)-0.5*std::log(2*dlib::pi)-log_sig_1c;

   // if(s0<2)  std::cout<<"log_ana_integrated_bi0 = "<<log_ana_integrated_bi0<<"\n";

  //std::cout<<"s0 = "<<s0<<"; log_prod_f_uij_given_bi0 = "<<log_prod_f_uij_given_bi0<<"\n";

   weighted_sum = weighted_sum + std::exp(log_prod_f_uij_given_bi0 + log_ana_integrated_bi0) * (*p_GH_points)(s0,1);

    }

       loglike = loglike + std::log(weighted_sum) -0.5*std::log(dlib::pi);
    }
 return (-1.0)*loglike;
 }

  // CHi-bar
  dbmat* p_cov_mle_inv; dbmat* p_mvn_sample_for_chi_bar_cdf;
  dbcolvec* p_one_mvn_sample_for_chi_bar_cdf;

   double obj_to_min_chi_bar_d4(const dbcolvec&theta)
  {
    return (trans(*p_one_mvn_sample_for_chi_bar_cdf-theta)*(*p_cov_mle_inv)*(*p_one_mvn_sample_for_chi_bar_cdf-theta) )(0,0);
  }

   dbcolvec grad_obj_to_min_chi_bar_d4(const dbcolvec&theta)
   {
     return (2*(*p_cov_mle_inv)*(*p_one_mvn_sample_for_chi_bar_cdf-theta) );
   }

  // double_num_chi_bar_sam
   void chi_bar_d4_cdf_weights(dbmat&cov_mle, dbmat&cov_mle_inv, dlib::rand&rd, int& num_chi_bar_sam, dbcolvec&proj_of_mle, dbmat&L_cov_mle, dbcolvec&std_mvn_sam,dbcolvec& L_mvn_sam, dbcolvec& lb_chi_bar_cdf_simu,dbcolvec& ub_chi_bar_cdf_simu,dbcolvec& weights)
   {
     double double_num_chi_bar_sam = double(num_chi_bar_sam);
     p_cov_mle_inv = &cov_mle_inv;
    //std::cout<<"*p_cov_mle_inv \n"<<*p_cov_mle_inv;
    //std::cout<<"\nL_cov_mle \n"<<L_cov_mle;
     int bfgs_fail = 0, s, k;
     int num_0_pos=0, num_1_pos=0, num_2_pos=0, num_3_pos=0, num_4_pos=0;
     int comp0_pos, comp1_pos, comp2_pos, comp3_pos, total_pos_comp;

     for(k=0; k<num_chi_bar_sam; ++k)
     {
         for(s=0;s<4;++s) { std_mvn_sam(s) = rd.get_random_gaussian(); }
         L_mvn_sam = L_cov_mle*std_mvn_sam;
         p_one_mvn_sample_for_chi_bar_cdf = &L_mvn_sam;

    proj_of_mle = 0; bfgs_fail = 0;
  find_min_box_constrainedConvFail(bfgs_search_strategy(),\
            objective_delta_stop_strategy(1e-5), \
            obj_to_min_chi_bar_d4, grad_obj_to_min_chi_bar_d4, proj_of_mle, lb_chi_bar_cdf_simu,ub_chi_bar_cdf_simu,bfgs_fail,2000);


    if(bfgs_fail==1) {std::cout<<"problem in chi-bar_cdf calculation";}
       else
       {
          comp0_pos = 0; comp1_pos = 0;  comp2_pos = 0; comp3_pos = 0;
          if(proj_of_mle(0) > 0 ) comp0_pos = 1;
          if(proj_of_mle(1) > 0 ) comp1_pos = 1;
          if(proj_of_mle(2) > 0 ) comp2_pos = 1;
          if(proj_of_mle(3) > 0 ) comp3_pos = 1;

      total_pos_comp = comp0_pos + comp1_pos + comp2_pos + comp3_pos;
       switch (total_pos_comp)
        {
        case 0:
           num_0_pos++;
            break;

          case 1:
          num_1_pos++;
            break;

         case 2:
          num_2_pos++;
            break;

         case 3:
          num_3_pos++;
            break;

           case 4:
        num_4_pos++;
            break;
       }

     }
     }
    weights(0) = double(num_0_pos)/double_num_chi_bar_sam;
    weights(1) = double(num_1_pos)/double_num_chi_bar_sam;
    weights(2) = double(num_2_pos)/double_num_chi_bar_sam;
    weights(3) = double(num_3_pos)/double_num_chi_bar_sam;
    weights(4) = double(num_4_pos)/double_num_chi_bar_sam;
   }

   double chi_bar_d4_cdf(const double&x, dbcolvec& weights)
   {
       double cdf = weights(0) + weights(1)*scythe::pchisq(x,1) + weights(2)*scythe::pchisq(x,2)+ weights(3)*scythe::pchisq(x,3) + weights(4)*scythe::pchisq(x,4);
       return cdf;
   }
 /*
// Laplace approximation for LME (exact)
    double g_bi_lme_d1(const dbcolvec& bi)
  {
   double loglike = 0, tij, residual, HTi, LTi;

  for(int j=*p_start_j_sub_i; j<(*p_endj_sub_i); ++j) //
  {
    tij = (*pDat_)(j,3);
    HTi = (*pDat_)(j,4);
    LTi = (*pDat_)(j,5);

            if (std::abs((*pDat_)(j,1)) <0.1)
             {

             }
             else
             {
              residual = (*pDat_)(j,2) - bi(0) - (*p_para_vec)(0)  - tij*((*p_para_vec)(1)) - ((*p_para_vec)(2)) * HTi - ((*p_para_vec)(3)) * LTi;
             loglike = loglike  -(*p_para_vec)(5) - 0.5*residual*residual/(std::exp(2*((*p_para_vec)(5)) ));
          }

   //Debug// DebugOut<<"j="<<j<<"; linear_pred_censored="<<linear_pred_censored<<"; residual="<<residual<<"; 0.5*residual*residual/(std::exp(2*((*p_para_vec)(8)) ))="<<0.5*residual*residual/(std::exp(2*((*p_para_vec)(8)) ))<<"; 0.5*(trans(ci)*(*p_cov_ranf_inv)*ci)(0,0)="<<0.5*(trans(ci)*(*p_cov_ranf_inv)*ci)(0,0)<<"\n";

  }
   // add density from random effects,  (*p_half_dim_rf)*std::log(2*dlib::pi)
   loglike = loglike - 0.5*std::log(2*dlib::pi) + 0.5*std::log(*p_det_cov_ranf_inv) -0.5*(trans(bi)*(*p_cov_ranf_inv)*bi)(0,0);
  return (-1.0)*loglike;
}

dbcolvec  Grad_g_bi_lme_d1(const dbcolvec& bi)
{
  int dim_ranef = bi.nr();
  dbcolvec Z_ijA_star(dim_ranef);
  dbcolvec  grad(dim_ranef);
  double tij, residual, HTi, LTi;

  grad = 0;
  for(int j=*p_start_j_sub_i; j<(*p_endj_sub_i); ++j) //
  {
    tij = (*pDat_)(j,3);
    HTi = (*pDat_)(j,4);
    LTi = (*pDat_)(j,5);
    Z_ijA_star(0) = 1;

          if (std::abs((*pDat_)(j,1)) <0.1)
             {

             }
             else
             {

              residual = (*pDat_)(j,2) - bi(0) - (*p_para_vec)(0)  - tij*((*p_para_vec)(1)) - ((*p_para_vec)(2)) * HTi - ((*p_para_vec)(3)) * LTi;
           grad = grad + (residual/(std::exp(2*((*p_para_vec)(5)) )) )*Z_ijA_star;
          }

  }
   // add density from random effects,
   grad = grad - (*p_cov_ranf_inv)*bi;

   return (-1.0)*grad;
}

dbmat  Hess_g_bi_lme_d1(const dbcolvec& bi)
{
  int dim_ranef = bi.nr();
  dbmat hess(dim_ranef,dim_ranef);
  dbcolvec Z_ijA_star(dim_ranef);
  double tij, residual, HTi, LTi;
  hess = 0;

   for(int j=*p_start_j_sub_i; j<(*p_endj_sub_i); ++j) //
  {
    tij = (*pDat_)(j,3);
    HTi = (*pDat_)(j,4);
    LTi = (*pDat_)(j,5);
    Z_ijA_star(0) = 1;


            if (std::abs((*pDat_)(j,1)) <0.1)
             {

             }
             else
             {

          residual = (*pDat_)(j,2) - bi(0) - (*p_para_vec)(0)  - tij*((*p_para_vec)(1)) - ((*p_para_vec)(2)) * HTi - ((*p_para_vec)(3)) * LTi;
          hess = hess - (1/(std::exp(2*((*p_para_vec)(5)) )) )*Z_ijA_star*trans(Z_ijA_star);

          }

  }
    hess = hess - (*p_cov_ranf_inv);

 return (-1.0)*hess;
}



  int Lap_neg_loglike_ctr_lme_d1 = 0;
  double Lap_neg_loglike_lme_d1(const dbcolvec& para_vec)
  {
     double loglike, loglike_i;
     int start_j_sub_i, ni, endj;

     int dim_rf = 1; double half_dim_rf = double(dim_rf)/2.0;
     dbcolvec cov_ranef_inv_Lvec(1);
     dbmat cov_ranef_inv_L(dim_rf,dim_rf), cov_ranef_inv(dim_rf,dim_rf);
     cov_ranef_inv_Lvec = subm(para_vec,range(4,4),range(0,0));
     EtaToL(cov_ranef_inv_Lvec,dim_rf,cov_ranef_inv_L);
    cov_ranef_inv =  cov_ranef_inv_L*trans(cov_ranef_inv_L);
    //Debug//  std::cout<<"In Lap, cov_ranef_inv_L\n"<<cov_ranef_inv_L;
  //Debug// DebugOut<<"cov_ranef_inv\n"<<cov_ranef_inv;



     double det_cov_ranf_inv = det(cov_ranef_inv_L)* det(cov_ranef_inv_L); //it is easier to compute det of a triangular matrix


    //Debug// std::cout<<"inside AGH neg loglike, cov_ranef_L (choleskey factor of cov_ranef):\n"<<cov_ranef_L<<"\n";

     dbcolvec bi_mode(dim_rf);
     dbmat SigHat(dim_rf,dim_rf); dbcolvec bi_mode_ub(dim_rf), bi_mode_lb(dim_rf), grad_ci(dim_rf);

     bi_mode_lb = -1.0*BIGNUM;
     bi_mode_ub = BIGNUM;


    int g_bi_mode_fail,  max_fun_eval = 200000; // max_iterations = 200,
    p_det_cov_ranf_inv = &det_cov_ranf_inv;
    p_cov_ranf_inv = &cov_ranef_inv;
    p_para_vec = &para_vec;

     loglike = 0;
     start_j_sub_i = 0;

    for(int i = 0; i < *p_nsub; ++i)  //Debug//  int i = 0; i<1; *p_nsub;
    {
         if(i > 0)
        { start_j_sub_i = start_j_sub_i \
         + int( (*p_ni_vec)(i-1) );  }

         ni = int( (*p_ni_vec)(i) );
         endj = start_j_sub_i + ni;

         p_start_j_sub_i = &start_j_sub_i;
         p_endj_sub_i = &endj;




   // compute ci_mode, detL
        g_bi_mode_fail = 0; bi_mode = 0.0;
  DebugOut<<"\n\n\ni="<<i<<"; numerical grad at bi_mode=0: "<<trans(gradcdif(g_bi_lme_d1,bi_mode))<<"i="<<i<<"; analytical grad at bi_mode=0: "<<trans(Grad_g_bi_lme_d1(bi_mode));
     find_minConvFail(bfgs_search_strategy(), objective_delta_stop_strategy(1e-7), g_bi_lme_d1, Grad_g_bi_lme_d1, bi_mode, -BIGNUM,g_bi_mode_fail,1000000);
   // Debug// .be_verbose(); std::cout<<"g_ci_mode_fail="<<g_ci_mode_fail<<"\n";
  //
  DebugOut<<"\ni="<<i<<"; numerical grad at final bi_mode: "<<trans(gradcdif(g_bi_lme_d1,bi_mode))<<"i="<<i<<"; analytical grad at final bi_mode: "<<trans(Grad_g_bi_lme_d1(bi_mode))<<"i="<<i<<"; numerical hess at final bi_mode: "<<HESScdif(g_bi_lme_d1,bi_mode)<<"i="<<i<<"; analytical hess at final bi_mode: "<<Hess_g_bi_lme_d1(bi_mode);

    // find_minConvFail(newton_search_strategy(Hess_g_jm_d5), objective_delta_stop_strategy(5e-7), g_jm_d5, Grad_g_jm_d5, ci_mode, -BIGNUM,g_ci_mode_fail,1000000);

    // find_min_trust_regionConvFail(objective_delta_stop_strategy(1e-11), cls_g_jm_d5(), ci_mode, g_ci_mode_fail, 300, 2);
   // find_min_bobyqaConvgFail(g_jm_d5,  ci_mode, 2*dim_rf+1,   ci_mode_lb, ci_mode_ub,  1.0,  1e-15,   max_fun_eval, g_ci_mode_fail);
   //

  //Debug// DebugOut<<"i="<<i<<"; ci_mode = "<<trans(ci_mode);
  //Debug// DebugOut<<"i="<<i<<"; numerical grad at final ci_mod: \n"<<trans(gradcdif(g_jm_d5,ci_mode,p_auto_step_size));
  //Debug// DebugOut<<"i="<<i<<"; analytical grad at final ci_mode: \n"<<trans(Grad_g_jm_d5(ci_mode));
    //Debug// DebugOut<<"i="<<i<<"; numerical hess at final ci_mode: \n"<<HESScdif(g_jm_d5,ci_mode,p_auto_step_size);
  //Debug// DebugOut<<"i="<<i<<"; analytical hess at final ci_mode: \n"<<Hess_g_jm_d5(ci_mode);

  //Debug//  DebugOut<<"i="<<i<<"; approximate analytical hess at final ci_mode: \n"<<ApHess_g_jm_d5(ci_mode);

      if(g_bi_mode_fail==1) {

  //  g_bi_mode_fail = 0; bi_mode = 0.0;
  // find_min_bobyqaConvgFail(g_bi,  bi_mode, 2*dim_rf+1,   bi_mode_lb, bi_mode_ub,  1.0,  1e-9,   max_fun_eval, g_bi_mode_fail);
    // std::cout<<"g_ci_mode_fail="<<g_ci_mode_fail<<"\n";

         if(g_bi_mode_fail==1)  {

     std::cout<<"Failure in the calculation of bi_mode \n";
  //Debug//   DebugOut<<std::setprecision(20)<<"Lap_neg_loglike_jm_d5_eval="<<Lap_neg_loglike_jm_d5_eval_noexp<<"; Failure in the calculation of ci_mode\n";
    return -BIGNUM;
       }
      }



      else {
        //Debug// SigHat = chol(HESScdif(g_jm_d5,ci_mode,p_auto_step_size));
             SigHat =   Hess_g_bi_lme_d1(bi_mode)  ;



    loglike_i =  -1.0*g_bi_lme_d1(bi_mode) + half_dim_rf*std::log(2*dlib::pi) + 0.5*std::log(det(SigHat));
    // std::cout<<"i="<<i<<"; g_jm_d5(ci_mode)="<<g_jm_d5(ci_mode)<<"; det(SigHat)="<<det(SigHat)<<"\n";
    // std::cout<<"i="<<i<<"; loglike_i="<<loglike_i<<"\n";
       loglike = loglike + loglike_i;

       }

    }

    Lap_neg_loglike_ctr_lme_d1++;
    //Debug// DebugOut<<"#like_eva="<<Lap_neg_loglike_jm_d5_eval_noexp<<"; loglike="<<loglike<<"\n";

   // if( (Lap_neg_loglike_jm_d5_eval_noexp%6000==0) ) {
 //std::cout<<"#like_eva="<<Lap_neg_loglike_jm_d5_eval_noexp<<"; neg_loglike="<<(-1.0)*loglike<<"\n";// <<"; para_vec: ";<<trans(para_vec)<<"\n"; // <<std::setprecision(20)   <<"; para_vec: "<<trans(para_vec)
  // }


 return (-1.0)*loglike;
  }


 dbcolvec Grad_Lap_neg_loglike_lme_d1(const dbcolvec& para_vec)
 {
     return gradcdif(Lap_neg_loglike_lme_d1,para_vec);
 }




 // parameterization is based on LOG-choleskey factorization of INVERSE covariance matrix
// requires pDat_, p_para_vec,  p_start_j_sub_i, p_endj_sub_i, p_cov_ranf_inv, and p_det_cov_ranf_inv to be READY
 double g_bi(const dbcolvec& bi)
  {
   double g = 0, tij, linear_pred_prob_posi, residual, HTi, LTi;

  for(int j=*p_start_j_sub_i; j<(*p_endj_sub_i); ++j) //
  {
    tij = (*pDat_)(j,3);
    HTi = (*pDat_)(j,4);
    LTi = (*pDat_)(j,5);
    linear_pred_prob_posi = bi(0) + (*p_para_vec)(0)  + tij*((*p_para_vec)(1)) + ((*p_para_vec)(2)) * HTi + ((*p_para_vec)(3)) * LTi ;
            if (std::abs((*pDat_)(j,1)) <0.1)
             {
                g = g - std::log(1 + std::exp(linear_pred_prob_posi));
             }
             else
             {
                 g = g + linear_pred_prob_posi- std::log(1 + std::exp(linear_pred_prob_posi));
              residual = (*pDat_)(j,2) - bi(1) - (*p_para_vec)(4)  - tij*((*p_para_vec)(5)) - ((*p_para_vec)(6)) * HTi - ((*p_para_vec)(7)) * LTi;
             g = g  -(*p_para_vec)(8) - 0.5*residual*residual/(std::exp(2*((*p_para_vec)(8)) ));
          }

   //Debug// DebugOut<<"j="<<j<<"; linear_pred_censored="<<linear_pred_censored<<"; residual="<<residual<<"; 0.5*residual*residual/(std::exp(2*((*p_para_vec)(8)) ))="<<0.5*residual*residual/(std::exp(2*((*p_para_vec)(8)) ))<<"; 0.5*(trans(ci)*(*p_cov_ranf_inv)*ci)(0,0)="<<0.5*(trans(ci)*(*p_cov_ranf_inv)*ci)(0,0)<<"\n";

  }
   // add density from random effects,  (*p_half_dim_rf)*std::log(2*dlib::pi)
   g = g - std::log(2*dlib::pi) + 0.5*std::log(*p_det_cov_ranf_inv) -0.5*(trans(bi)*(*p_cov_ranf_inv)*bi)(0,0);
  return (-1.0)*g;
}

dbcolvec  Grad_g_bi(const dbcolvec& bi)
{
  int dim_ranef = bi.nr();
  dbcolvec Z_ijA(dim_ranef), Z_ijA_star(dim_ranef);
  dbcolvec  grad(dim_ranef);
  double tij, linear_pred_prob_posi, residual, HTi, LTi;

  grad = 0;
  for(int j=*p_start_j_sub_i; j<(*p_endj_sub_i); ++j) //
  {
    tij = (*pDat_)(j,3);
    HTi = (*pDat_)(j,4);
    LTi = (*pDat_)(j,5);
    Z_ijA = 0; Z_ijA_star = 0;
    Z_ijA(0) = 1; Z_ijA_star(1) = 1;

    linear_pred_prob_posi = bi(0) + (*p_para_vec)(0)  + tij*((*p_para_vec)(1)) + ((*p_para_vec)(2)) * HTi + ((*p_para_vec)(3)) * LTi;

          if (std::abs((*pDat_)(j,1)) <0.1)
             {
             grad = grad -1.0*( 1/( 1+std::exp(-1.0*linear_pred_prob_posi ) )  )*Z_ijA;
             }
             else
             {
           grad = grad + (1- 1/( 1+std::exp(-1.0*linear_pred_prob_posi ) )  )*Z_ijA;
              residual = (*pDat_)(j,2) - bi(1) - (*p_para_vec)(4)  - tij*((*p_para_vec)(5)) - ((*p_para_vec)(6)) * HTi - ((*p_para_vec)(7)) * LTi;
           grad = grad + (residual/(std::exp(2*((*p_para_vec)(8)) )) )*Z_ijA_star;
          }

  }
   // add density from random effects,
   grad = grad - (*p_cov_ranf_inv)*bi;

   return (-1.0)*grad;
}

dbmat  Hess_g_bi(const dbcolvec& bi)
{
  int dim_ranef = bi.nr();
  dbmat hess(dim_ranef,dim_ranef);
  dbcolvec Z_ijA(dim_ranef), Z_ijA_star(dim_ranef);
  double tij, linear_pred_prob_posi, residual, HTi, LTi, temp;
  hess = 0;

   for(int j=*p_start_j_sub_i; j<(*p_endj_sub_i); ++j) //
  {
    tij = (*pDat_)(j,3);
    HTi = (*pDat_)(j,4);
    LTi = (*pDat_)(j,5);
    Z_ijA = 0; Z_ijA_star = 0;
    Z_ijA(0) = 1; Z_ijA_star(1) = 1;

    linear_pred_prob_posi = bi(0) + (*p_para_vec)(0)  + tij*((*p_para_vec)(1)) + ((*p_para_vec)(2)) * HTi + ((*p_para_vec)(3)) * LTi ;
    temp = 1 + std::exp(-1.0*linear_pred_prob_posi);

            if (std::abs((*pDat_)(j,1)) <0.1)
             {
                hess = hess -1.0*((temp-1)/(temp*temp))*Z_ijA*trans(Z_ijA);
             }
             else
             {
          hess = hess -1.0*((temp-1)/(temp*temp))*Z_ijA*trans(Z_ijA);
          residual = (*pDat_)(j,2) - bi(1) - (*p_para_vec)(4)  - tij*((*p_para_vec)(5)) - ((*p_para_vec)(6)) * HTi - ((*p_para_vec)(7)) * LTi;
          hess = hess - (1/(std::exp(2*((*p_para_vec)(8)) )) )*Z_ijA_star*trans(Z_ijA_star);

          }

  }
    hess = hess - (*p_cov_ranf_inv);

 return (-1.0)*hess;
}



  int Lap_neg_loglike_ctr = 0;
  double Lap_neg_loglike(const dbcolvec& para_vec)
  {
     double loglike, loglike_i;
     int start_j_sub_i, ni, endj;

     int dim_rf = 2; double half_dim_rf = double(dim_rf)/2.0;
     dbcolvec cov_ranef_inv_Lvec(3);
     dbmat cov_ranef_inv_L(2,2), cov_ranef_inv(2,2);
     cov_ranef_inv_Lvec = subm(para_vec,range(9,11),range(0,0));
     EtaToL(cov_ranef_inv_Lvec,dim_rf,cov_ranef_inv_L);
    cov_ranef_inv =  cov_ranef_inv_L*trans(cov_ranef_inv_L);
    //Debug//  std::cout<<"In Lap, cov_ranef_inv_L\n"<<cov_ranef_inv_L;
  //Debug// DebugOut<<"cov_ranef_inv\n"<<cov_ranef_inv;



     double det_cov_ranf_inv = det(cov_ranef_inv_L)* det(cov_ranef_inv_L); //it is easier to compute det of a triangular matrix


    //Debug// std::cout<<"inside AGH neg loglike, cov_ranef_L (choleskey factor of cov_ranef):\n"<<cov_ranef_L<<"\n";

     dbcolvec bi_mode(dim_rf);
     dbmat SigHat(dim_rf,dim_rf); dbcolvec bi_mode_ub(dim_rf), bi_mode_lb(dim_rf), grad_ci(dim_rf);

     bi_mode_lb = -1.0*BIGNUM;
     bi_mode_ub = BIGNUM;


    int g_bi_mode_fail,  max_fun_eval = 200000; // max_iterations = 200,
    p_det_cov_ranf_inv = &det_cov_ranf_inv;
    p_cov_ranf_inv = &cov_ranef_inv;
    p_para_vec = &para_vec;

     loglike = 0;
     start_j_sub_i = 0;

    for(int i = 0; i < *p_nsub; ++i)  //Debug//  int i = 0; i<1; *p_nsub;
    {
         if(i > 0)
        { start_j_sub_i = start_j_sub_i \
         + int( (*p_ni_vec)(i-1) );  }

         ni = int( (*p_ni_vec)(i) );
         endj = start_j_sub_i + ni;

         p_start_j_sub_i = &start_j_sub_i;
         p_endj_sub_i = &endj;




   // compute ci_mode, detL
        g_bi_mode_fail = 0; bi_mode = 0.0;
  DebugOut<<"i="<<i<<"; numerical grad at bi_mode=0: "<<trans(gradcdif(g_bi,bi_mode))<<"i="<<i<<"; analytical grad at bi_mode=0: "<<trans(Grad_g_bi(bi_mode));
     find_minConvFail(bfgs_search_strategy(), objective_delta_stop_strategy(1e-7), g_bi, Grad_g_bi, bi_mode, -BIGNUM,g_bi_mode_fail,1000000);
   // Debug// .be_verbose(); std::cout<<"g_ci_mode_fail="<<g_ci_mode_fail<<"\n";
  //  DebugOut<<"i="<<i<<"; numerical grad at final bi_mode: "<<trans(gradcdif(g_bi,bi_mode))<<"i="<<i<<"; analytical grad at final bi_mode: "<<trans(Grad_g_bi(bi_mode))<<"i="<<i<<"; numerical hess at final bi_mode: "<<HESScdif(g_bi,bi_mode)<<"i="<<i<<"; analytical hess at final bi_mode: "<<Hess_g_bi(bi_mode);

    // find_minConvFail(newton_search_strategy(Hess_g_jm_d5), objective_delta_stop_strategy(5e-7), g_jm_d5, Grad_g_jm_d5, ci_mode, -BIGNUM,g_ci_mode_fail,1000000);

    // find_min_trust_regionConvFail(objective_delta_stop_strategy(1e-11), cls_g_jm_d5(), ci_mode, g_ci_mode_fail, 300, 2);
   // find_min_bobyqaConvgFail(g_jm_d5,  ci_mode, 2*dim_rf+1,   ci_mode_lb, ci_mode_ub,  1.0,  1e-15,   max_fun_eval, g_ci_mode_fail);
   //

  //Debug// DebugOut<<"i="<<i<<"; ci_mode = "<<trans(ci_mode);
  //Debug// DebugOut<<"i="<<i<<"; numerical grad at final ci_mod: \n"<<trans(gradcdif(g_jm_d5,ci_mode,p_auto_step_size));
  //Debug// DebugOut<<"i="<<i<<"; analytical grad at final ci_mode: \n"<<trans(Grad_g_jm_d5(ci_mode));
    //Debug// DebugOut<<"i="<<i<<"; numerical hess at final ci_mode: \n"<<HESScdif(g_jm_d5,ci_mode,p_auto_step_size);
  //Debug// DebugOut<<"i="<<i<<"; analytical hess at final ci_mode: \n"<<Hess_g_jm_d5(ci_mode);

  //Debug//  DebugOut<<"i="<<i<<"; approximate analytical hess at final ci_mode: \n"<<ApHess_g_jm_d5(ci_mode);

      if(g_bi_mode_fail==1) {

    g_bi_mode_fail = 0; bi_mode = 0.0;
  find_min_bobyqaConvgFail(g_bi,  bi_mode, 2*dim_rf+1,   bi_mode_lb, bi_mode_ub,  1.0,  1e-9,   max_fun_eval, g_bi_mode_fail);
    // std::cout<<"g_ci_mode_fail="<<g_ci_mode_fail<<"\n";

         if(g_bi_mode_fail==1)  {

     std::cout<<"Failure in the calculation of bi_mode \n";
  //Debug//   DebugOut<<std::setprecision(20)<<"Lap_neg_loglike_jm_d5_eval="<<Lap_neg_loglike_jm_d5_eval_noexp<<"; Failure in the calculation of ci_mode\n";
    return -BIGNUM;
       }
      }



      else {
        //Debug// SigHat = chol(HESScdif(g_jm_d5,ci_mode,p_auto_step_size));
             SigHat =   Hess_g_bi(bi_mode)  ;



    loglike_i =  -1.0*g_bi(bi_mode) + half_dim_rf*std::log(2*dlib::pi) - 0.5*std::log(det(SigHat));
    // std::cout<<"i="<<i<<"; g_jm_d5(ci_mode)="<<g_jm_d5(ci_mode)<<"; det(SigHat)="<<det(SigHat)<<"\n";
    // std::cout<<"i="<<i<<"; loglike_i="<<loglike_i<<"\n";
       loglike = loglike + loglike_i;

       }

    }

    Lap_neg_loglike_ctr++;
    //Debug// DebugOut<<"#like_eva="<<Lap_neg_loglike_jm_d5_eval_noexp<<"; loglike="<<loglike<<"\n";

   // if( (Lap_neg_loglike_jm_d5_eval_noexp%6000==0) ) {
 //std::cout<<"#like_eva="<<Lap_neg_loglike_jm_d5_eval_noexp<<"; neg_loglike="<<(-1.0)*loglike<<"\n";// <<"; para_vec: ";<<trans(para_vec)<<"\n"; // <<std::setprecision(20)   <<"; para_vec: "<<trans(para_vec)
  // }


 return (-1.0)*loglike;
  }


 dbcolvec Grad_Lap_neg_loglike(const dbcolvec& para_vec)
 {
     return gradcdif(Lap_neg_loglike,para_vec);
 }

  */







#endif // HEADER_H_INCLUDED
