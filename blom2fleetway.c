/*
**  Transform Blom survey coordinates to Fleetway.
*/
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

double utran[3][3]={1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0}, rtran[3][3];
#define rotation_angle 0.371*M_PI/180.0
#define x_offset 0.399
#define y_offset -0.523

static void
vsprd(double *a1, double *a2, int m, int n, double *a3)
{
        int	i,j,tn;
        double  a1s, *at1, *at2, *at3, *att2;

        if (n <= 0 || m <= 0)
            return;
        for (i = m; i>0; i--) {
        /* scaler-vector multiplication */
            for (tn=n,a1s= *a1,at2=a2,at3=a3;tn>0;tn--) {
                *at3    = a1s * (*at2);
                at2     += m; at3 += m;
            }
            for (at1=a1+m,at2=a2+1,j=2; j<=m; j++) {
            /* vpiv  -- vector pivot operation */
                for (a1s= *at1,tn=n,att2=at2,at3=a3;tn>0;tn--) {
                    *at3        += a1s * (*att2);
                    att2        += m; at3 += m;
                }
                at1     += m; at2++;
            }
            a1++; a3++;
        }
        return;
}

/* scalar-vector multiplication */
static void
vsmy(double s, double *a1, int incr1, double *a2, int incr2, int n)
{

        for (; n>0; n--) {
            *a2 = s * *a1;
            a1  += incr1;
            a2  += incr2;
        }
}

static void
init(void)
{
	double xform[3][3], ttran[3][3], sinx, cosx, f;

    /*
    **  Setup the rotation transform.
    */
	memset(xform, 0, sizeof (xform));
	xform[2][2] = 1.0;
	sinx = sin(rotation_angle);
	cosx = cos(rotation_angle);
	xform[0][0] = cosx;
	xform[0][1] = sinx;
	xform[1][0] = -sinx;
	xform[1][1] = cosx;
	vsprd(xform, utran, 3, 3, ttran);
	memmove(utran, ttran, sizeof(ttran));
    /*
    **  Setup the translation transform.
    */
	memset(xform, 0, sizeof (xform));
	xform[0][0] = xform[1][1] = xform[2][2] =  1.0;
	xform[2][0] = x_offset;
	xform[2][1] = y_offset;
	vsprd(xform, utran, 3, 3, ttran);
	memmove(utran, ttran, sizeof(ttran));
    /*
    **  Compute the inverse transform.
    */
	memmove(rtran, ttran, sizeof(rtran));
	rtran[0][0] =  ttran[1][1];
	rtran[1][0] = -ttran[1][0];
	rtran[2][0] =  ttran[2][1]*ttran[1][0]-ttran[1][1]*ttran[2][0];
	rtran[0][1] = -ttran[0][1];
	rtran[1][1] =  ttran[0][0];
	rtran[2][1] =  ttran[2][0]*ttran[0][1]-ttran[2][1]*ttran[0][0];

	f           =  ttran[0][0]*ttran[1][1]-ttran[0][1]*ttran[1][0];
	f           =  1.0/f;
	vsmy(f,rtran,3,rtran,3,3);
	/* scaler-vector multiplication */
	vsmy(f,rtran+1,3,rtran+1,3,3);
	/* scaler-vector multiplication */
}

static void
b2f(double xv, double yv, double *xa, double *ya)
{
	*xa = xv*utran[0][0] + yv*utran[1][0] + utran[2][0];
	*ya = xv*utran[0][1] + yv*utran[1][1] + utran[2][1];
}

static void
f2b(double xv, double yv, double *xa, double *ya)
{
	*xa = xv*rtran[0][0] + yv*rtran[1][0] + rtran[2][0];
	*ya = xv*rtran[0][1] + yv*rtran[1][1] + rtran[2][1];
}

int
main(int argc, char **argv)
{
	double	xv, yv, xa, ya;
	char	line[80];
	int	b2fcvt;

	if (argc != 2 || (strcmp(argv[1],"b2f") && strcmp(argv[1], "f2b"))) {
	    fprintf(stderr, "Usage: b2f b2f|f2b\n");
	    exit(1);
	}
	if (!strcmp(argv[1], "b2f"))
	    b2fcvt = 1;
	else
	    b2fcvt = 0;
	init();

	while (1) {
	    if (fgets(line, sizeof (line), stdin) == NULL)
		break;
	    if (sscanf(line, "%lf %lf", &xv, &yv) != 2)
		break;
	    if (b2fcvt)
		b2f(xv, yv, &xa, &ya);
	    else
		f2b(xv, yv, &xa, &ya);
	    printf("%.3f %.3f\n", xa, ya);
	}
	return (0);
}
