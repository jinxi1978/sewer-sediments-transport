//-----------------------------------------------------------------------------
//   xsect.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//             03/14/17   (Build 5.1.012)
//             11/27/17   (Build 5.1.013)
//   Author:   L. Rossman (EPA)
//             M. Tryby (EPA)
//
//   Cross section geometry functions.
//
//   The primary functions are:
//      getAofY   -- returns area given depth
//      getWofY   -- returns top width given depth
//      getRofY   -- returns hyd. radius given depth
//      getYofA   -- returns flow depth given area
//      getRofA   -- returns hyd. radius given area
//      getSofA   -- returns section factor given area
//      getAofS   -- returns area given section factor
//      getdSdA   -- returns derivative of section factor w.r.t. area
//   where
//      Y = flow depth
//      A = flow area
//      R = hyd. radius
//      S = section factor = A*R^(2/3)
//
//   Build 5.1.012:
//   - Height at max. width for Modified Baskethandle shape corrected.
//
//   Build 5.1.013:
//   - Width at full height set to 0 for closed rectangular shape.
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <math.h>
#include "headers.h"
#include "findroot.h"

#include "xsect.h"    // File containing geometry tables for rounded shapes
using namespace SWMMCPP;

//=============================================================================
double xsectClass::lookup(double x, double *table, int nItems)
//
//  Input:   x = value of independent variable in a geometry table
//           table = ptr. to geometry table
//           nItems = number of equally spaced items in table
//  Output:  returns value of dependent table variable
//  Purpose: looks up a value in a geometry table (i.e., finds y given x).
//
{
	double  delta, x0, x1, y, y2;
	int     i;

	// --- find which segment of table contains x
	delta = 1.0 / (nItems - 1);
	i = (int)(x / delta);
	if (i >= nItems - 1) return table[nItems - 1];

	// --- compute x at start and end of segment
	x0 = i * delta;
	x1 = (i + 1) * delta;

	// --- linearly interpolate a y-value
	y = table[i] + (x - x0) * (table[i + 1] - table[i]) / delta;

	// --- use quadratic interpolation for low x value
	if (i < 2)
	{
		y2 = y + (x - x0) * (x - x1) / (delta*delta) *
			(table[i] / 2.0 - table[i + 1] + table[i + 2] / 2.0);
		if (y2 > 0.0) y = y2;
	}
	if (y < 0.0) y = 0.0;
	return y;
}

double xsectClass::rect_closed_getRofA(double a)
{
	double p;
	if (a <= 0.0)   return 0.0;
	p = m_wMax + 2.*a / m_wMax; // Wetted Perim = width + 2*area/width
	if (a / m_aFull > RECT_ALFMAX)
	{
		p += (a / m_aFull - RECT_ALFMAX) / (1.0 - RECT_ALFMAX) * m_wMax;
	}
	return a / p;
}

double xsectClass::rect_triang_getYofA(double a)
{
	// below upper section
	if (a <= m_aBot) return sqrt(a / m_sBot);

	// above bottom section
	else return m_yBot + (a - m_aBot) / m_wMax;
}

double xsectClass::rect_triang_getRofA(double a)
{
	double y;
	double p, alf;

	if (a <= 0.0)   return 0.0;
	y = rect_triang_getYofA(a);

	// below upper section
	if (y <= m_yBot) return a / (2. * y * m_rBot);

	// wetted perimeter without contribution of top surface
	p = 2. * m_yBot * m_rBot + 2. * (y - m_yBot);

	// top-surface contribution
	alf = (a / m_aFull) - RECT_TRIANG_ALFMAX;
	if (alf > 0.0) p += alf / (1.0 - RECT_TRIANG_ALFMAX) * m_wMax;
	return a / p;
}

double xsectClass::getThetaOfAlpha(double alpha)
{
	int    k;
	double theta, theta1, ap, d;

	if (alpha > 0.04) theta = 1.2 + 5.08 * (alpha - 0.04) / 0.96;
	else theta = 0.031715 - 12.79384 * alpha + 8.28479 * sqrt(alpha);
	theta1 = theta;
	ap = (2.0*PI) * alpha;
	for (k = 1; k <= 40; k++)
	{
		d = -(ap - theta + sin(theta)) / (1.0 - cos(theta));
		// --- modification to improve convergence for large theta
		if (d > 1.0) d = SIGN(1.0, d);
		theta = theta - d;
		if (fabs(d) <= 0.0001) return theta;
	}
	return theta1;
}

double xsectClass::getYcircular(double alpha)
{
	double theta;
	if (alpha >= 1.0) return 1.0;
	if (alpha <= 0.0) return 0.0;
	if (alpha <= 1.0e-5)
	{
		theta = pow(37.6911*alpha, 1. / 3.);
		return theta * theta / 16.0;
	}
	theta = getThetaOfAlpha(alpha);
	return (1.0 - cos(theta / 2.)) / 2.0;
}

double xsectClass::rect_round_getYofA(double a)
{
	double alpha;

	// --- if above circular bottom:
	if (a > m_aBot)
		return m_yBot + (a - m_aBot) / m_wMax;

	// --- otherwise use circular xsection method to find height
	alpha = a / (PI * m_rBot * m_rBot);
	if (alpha < 0.04) return (2.0 * m_rBot) * getYcircular(alpha);
	return (2.0 * m_rBot) * lookup(alpha, Y_Circ, N_Y_Circ);
}

double xsectClass::rect_round_getRofA(double a)
{
	double y1, theta1, p, arg;

	// --- if above circular invert ...
	if (a <= 0.0) return 0.0;
	if (a > m_aBot)
	{
		// wetted perimeter without contribution of top surface
		y1 = (a - m_aBot) / m_wMax;
		theta1 = 2.0 * asin(m_wMax / 2.0 / m_rBot);
		p = m_rBot*theta1 + 2.0*y1;

		// top-surface contribution
		arg = (a / m_aFull) - RECT_ROUND_ALFMAX;
		if (arg > 0.0) p += arg / (1.0 - RECT_ROUND_ALFMAX) * m_wMax;
		return a / p;
	}

	// --- if within circular invert ...
	y1 = rect_round_getYofA(a);
	theta1 = 2.0*acos(1.0 - y1 / m_rBot);
	p = m_rBot * theta1;
	return a / p;
}

int xsectClass::xsect_isOpen()
//
//  Input:   type = type of xsection shape
//  Output:  returns 1 if xsection is open, 0 if not
//  Purpose: determines if a xsection type is open or closed.
//
{
    return ((Amax[m_type] >= 1.0) ? 1 : 0);
}

double xsectClass::getScircular(double alpha)
{
	double theta;
	if (alpha >= 1.0) return 1.0;
	if (alpha <= 0.0) return 0.0;
	if (alpha <= 1.0e-5)
	{
		theta = pow(37.6911*alpha, 1. / 3.);
		return pow(theta, 13. / 3.) / 124.4797;
	}
	theta = getThetaOfAlpha(alpha);
	return pow((theta - sin(theta)), 5. / 3.) / (2.0 * PI) / pow(theta, 2. / 3.);
}

double xsectClass::circ_getSofA(double a)
{
	double alpha = a / m_aFull;

	// --- use special function for small a/aFull
	if (alpha < 0.04) return m_sFull * getScircular(alpha);

	// --- otherwise use table
	else
		return m_sFull * lookup(alpha, S_Circ, N_S_Circ);
}

double xsectClass::rect_open_getSofA(double a)
{
	double y = a / m_wMax;
	double r = a / ((2.0 - m_sBot)*y + m_wMax);
	return a * pow(r, 2. / 3.);
}

double xsectClass::rect_triang_getSofA(double a)
{
	// --- if a > area corresponding to sMax, then
	//     interpolate between sMax and Sfull
	double alfMax = RECT_TRIANG_ALFMAX;
	if (a / m_aFull > alfMax)
		return m_sMax + (m_sFull - m_sMax) *
		(a / m_aFull - alfMax) / (1.0 - alfMax);

	// --- otherwise use regular formula
	else return a * pow(rect_triang_getRofA(a), 2. / 3.);
}

double xsectClass::rect_round_getSofA(double a)
{
	double alpha, aFull, sFull;

	// --- if a > area corresponding to sMax,
	//     interpolate between sMax and sFull
	double alfMax = RECT_ROUND_ALFMAX;
	if (a / m_aFull > alfMax)
	{
		return m_sMax + (m_sFull - m_sMax) *
			(a / m_aFull - alfMax) / (1.0 - alfMax);
	}

	// --- if above circular invert, use generic function
	else if (a > m_aBot)
	{
		return a * pow(xsect_getRofA(a), 2. / 3.);
	}

	// --- otherwise use circular xsection function applied
	//     to full circular shape of bottom section
	else
	{
		aFull = PI * m_rBot * m_rBot;
		alpha = a / aFull;
		sFull = m_sBot;

		// --- use special function for small a/aFull
		if (alpha < 0.04) return sFull * getScircular(alpha);

		// --- otherwise use table
		else return sFull * lookup(alpha, S_Circ, N_S_Circ);
	}
}

double xsectClass::xsect_getSofA(double a)
//
//  Input:   xsect = ptr. to a cross section data structure
//           a = area (ft2)
//  Output:  returns section factor (ft^(8/3))
//  Purpose: computes xsection's section factor at a given area.
//
{
	double alpha = a / m_aFull;
	double r;
	switch (m_type)
	{
	case FORCE_MAIN:
	case CIRCULAR:
		return circ_getSofA(a);

	case EGGSHAPED:
		return m_sFull * lookup(alpha, S_Egg, N_S_Egg);

	case HORSESHOE:
		return m_sFull * lookup(alpha, S_Horseshoe, N_S_Horseshoe);

	case GOTHIC:
		return m_sFull * lookup(alpha, S_Gothic, N_S_Gothic);

	case CATENARY:
		return m_sFull * lookup(alpha, S_Catenary, N_S_Catenary);

	case SEMIELLIPTICAL:
		return m_sFull * lookup(alpha, S_SemiEllip, N_S_SemiEllip);

	case BASKETHANDLE:
		return m_sFull * lookup(alpha, S_BasketHandle, N_S_BasketHandle);

	case SEMICIRCULAR:
		return m_sFull * lookup(alpha, S_SemiCirc, N_S_SemiCirc);

	case RECT_CLOSED:
		return rect_closed_getSofA(a);

	case RECT_OPEN:
		return rect_open_getSofA(a);

	case RECT_TRIANG:
		return rect_triang_getSofA(a);

	case RECT_ROUND:
		return rect_round_getSofA(a);

	default:
		if (a == 0.0) return 0.0;
		r = xsect_getRofA(a);
		if (r < TINY) return 0.0;
		return a * pow(r, 2. / 3.);
	}
}

double xsectClass::filled_circ_getRofY( double y)
{
	double a, r, p;

	// --- temporarily remove filled portion of circle
	m_yFull += m_yBot;
	m_aFull += m_aBot;
	y += m_yBot;

	// --- get area,  hyd. radius & wetted perimeter of unfilled circle
	a = circ_getAofY(y);
	r = 0.25 * m_yFull * lookup(y / m_yFull, R_Circ, N_R_Circ);
	p = (a / r);

	// --- reduce area and wetted perimeter by amount of filled circle
	//     (rBot = filled perimeter, sBot = filled width)
	a = a - m_aBot;
	p = p - m_rBot + m_sBot;

	// --- compute actual hyd. radius & restore xsect parameters
	r = a / p;
	m_yFull -= m_yBot;
	m_aFull -= m_aBot;
	return r;
}

double xsectClass::rect_triang_getRofY(double y)
{
	double p, a, alf;

	// y is below upper rectangular section
	if (y <= m_yBot) return y * m_sBot / (2. * m_rBot);

	// area
	a = m_aBot + (y - m_yBot) * m_wMax;

	// wetted perimeter without contribution of top surface
	p = 2. * m_yBot * m_rBot + 2. * (y - m_yBot);

	// top-surface contribution
	alf = (a / m_aFull) - RECT_TRIANG_ALFMAX;
	if (alf > 0.0) p += alf / (1.0 - RECT_TRIANG_ALFMAX) * m_wMax;
	return a / p;
}

double xsectClass::rect_round_getAofY(double y)
{
	double theta1;

	// --- if above circular invert...
	if (y > m_yBot)
		return m_aBot + (y - m_yBot) * m_wMax;

	// --- find area of circular section
	theta1 = 2.0*acos(1.0 - y / m_rBot);
	return 0.5 * m_rBot * m_rBot * (theta1 - sin(theta1));
}

double xsectClass::rect_round_getRofY(double y)
{
	double theta1;

	// --- if above top of circular bottom, use RofA formula
	if (y <= 0.0) return 0.0;
	if (y > m_yBot)
		return rect_round_getRofA(rect_round_getAofY(y));

	// --- find hyd. radius of circular section
	theta1 = 2.0*acos(1.0 - y / m_rBot);
	return 0.5 * m_rBot * (1.0 - sin(theta1)) / theta1;
}

double xsectClass::trapez_getAofY(double y)
{
	return (m_yBot + m_sBot * y) * y;
}

double xsectClass::trapez_getRofY(double y)
{
	if (y == 0.0) return 0.0;
	return trapez_getAofY(y) / (m_yBot + y * m_rBot);
}

double xsectClass::triang_getRofY(double y)
{
	return (y * m_sBot) / (2. * m_rBot);
}

double xsectClass::parab_getYofA(double a)
{
	return pow((3. / 4.) * a / m_rBot, 2. / 3.);
}

double xsectClass::parab_getRofA( double a)
{
	if (a <= 0.0) return 0.0;
	return a / parab_getPofY(parab_getYofA(a));
}

double xsectClass::parab_getPofY(double y)
{
	double x = 2. * sqrt(y) / m_rBot;
	double t = sqrt(1.0 + x * x);
	return 0.5 * m_rBot * m_rBot * (x * t + log(x + t));
}

double xsectClass::parab_getAofY(double y)
{
	return (4. / 3. * m_rBot * y * sqrt(y));
}

double xsectClass::parab_getRofY(double y)
{
	if (y <= 0.0) return 0.0;
	return parab_getAofY(y) / parab_getPofY(y);
}

double xsectClass::parab_getWofY(double y)
{
	return 2.0 * m_rBot * sqrt(y);
}

double xsectClass::powerfunc_getYofA(double a)
{
	return pow(a / m_rBot, 1.0 / (m_sBot + 1.0));
}

double xsectClass::powerfunc_getRofA(double a)
{
	if (a <= 0.0) return 0.0;
	return a / powerfunc_getPofY(powerfunc_getYofA(a));
}

double xsectClass::powerfunc_getPofY(double y)
{
	double dy1 = 0.02 * m_yFull;
	double h = (m_sBot + 1.0) * m_rBot / 2.0;
	double m = m_sBot;
	double p = 0.0;
	double y1 = 0.0;
	double x1 = 0.0;
	double x2, y2, dx, dy;
	do
	{
		y2 = y1 + dy1;
		if (y2 > y) y2 = y;
		x2 = h * pow(y2, m);
		dx = x2 - x1;
		dy = y2 - y1;
		p += sqrt(dx*dx + dy*dy);
		x1 = x2;
		y1 = y2;
	} while (y2 < y);
	return 2.0 * p;
}

double xsectClass::powerfunc_getAofY(double y)
{
	return m_rBot * pow(y, m_sBot + 1.0);
}

double xsectClass::powerfunc_getRofY(double y)
{
	if (y <= 0.0) return 0.0;
	return powerfunc_getAofY(y) / powerfunc_getPofY(y);
}

double xsectClass::powerfunc_getWofY(double y)
{
	return (m_sBot + 1.0) * m_rBot * pow(y, m_sBot);
}

double xsectClass::filled_circ_getAofY(double y)
{
	double a;

	// --- temporarily remove filled portion of circle
	m_yFull += m_yBot;
	m_aFull += m_aBot;
	y += m_yBot;

	// --- find area of unfilled circle
	a = circ_getAofY(y);

	// --- restore original values
	a -= m_aBot;
	m_yFull -= m_yBot;
	m_aFull -= m_aBot;
	return a;
}

int xsectClass::locate(double y, double *table, int jLast)
//
//  Input:   y      = value being located in table
//           table  = ptr. to table with monotonically increasing entries
//           jLast  = highest table entry index to search over
//  Output:  returns index j of table such that table[j] <= y <= table[j+1]
//  Purpose: uses bisection method to locate the highest table index whose
//           table entry does not exceed a given value.
//
//  Notes:   This function is only used in conjunction with invLookup().
//
{
	int j;
	int j1 = 0;
	int j2 = jLast;

	// Check if value <= first table entry
	if (y <= table[0]) return 0;

	// Check if value >= the last entry
	if (y >= table[jLast]) return jLast;

	// While a portion of the table still remains
	while (j2 - j1 > 1)
	{
		// Find midpoint of remaining portion of table
		j = (j1 + j2) >> 1;

		// Value is greater or equal to midpoint: search from midpoint to j2
		if (y >= table[j]) j1 = j;

		// Value is less than midpoint: search from j1 to midpoint
		else j2 = j;
	}

	// Return the lower index of the remaining interval,
	return j1;
}

double xsectClass::invLookup(double y, double *table, int nItems)
//
//  Input:   y = value of dependent variable in a geometry table
//           table = ptr. to geometry table
//           nItems = number of equally spaced items in table
//  Output:  returns value of independent table variable
//  Purpose: performs inverse lookup in a geometry table (i.e., finds
//           x given y).
//
//  Notes:   This function assumes that the geometry table has either strictly
//           increasing entries or that the maximum entry is always third
//           from the last (which is true for all section factor tables). In
//           the latter case, the location of a large y can be ambiguous
//           -- it can be both below and above the location of the maximum.
//           In such cases this routine searches only the interval above
//           the maximum (i.e., the last 2 segments of the table).
//
//           nItems-1 is the highest subscript for the table's data.
//
//           The x value's in a geometry table lie between 0 and 1.
//
{
	double dx;               // x-increment of table
	double x, x0, dy;        // interpolation variables
	int    n;                // # items in increasing portion of table
	int    i;                // lower table index that brackets y

							 // --- compute table's uniform x-increment
	dx = 1.0 / (double)(nItems - 1);

	// --- truncate item count if last 2 table entries are decreasing
	n = nItems;
	if (table[n - 3] > table[n - 1]) n = n - 2;

	// --- check if y falls in decreasing portion of table
	if (n < nItems && y > table[nItems - 1])
	{
		if (y >= table[nItems - 3]) return (n - 1) * dx;
		if (y <= table[nItems - 2]) i = nItems - 2;
		else i = nItems - 3;
	}

	// --- otherwise locate the interval where y falls in the table
	else i = locate(y, table, n - 1);
	if (i >= n - 1) return (n - 1) * dx;

	// --- compute x at start and end of segment
	x0 = i * dx;

	// --- linearly interpolate an x value
	dy = table[i + 1] - table[i];
	if (dy == 0.0) x = x0;
	else x = x0 + (y - table[i]) * dx / dy;
	if (x < 0.0) x = 0.0;
	if (x > 1.0) x = 1.0;
	return x;
}

double xsectClass::generic_getdSdA(double a)
//
//  Input:   xsect = ptr. to cross section data structure
//           a = area (ft2)
//  Output:  returns derivative of section factor w.r.t. area (ft^2/3)
//  Purpose: computes derivative of section factor w.r.t area
//           using central difference approximation.
//
{
	double a1, a2;
	double alpha = a / m_aFull;
	double alpha1 = alpha - 0.001;
	double alpha2 = alpha + 0.001;
	if (alpha1 < 0.0) alpha1 = 0.0;
	a1 = alpha1 * m_aFull;
	a2 = alpha2 * m_aFull;
	return (xsect_getSofA(a2) - xsect_getSofA(a1)) / (a2 - a1);
}

double xsectClass::rect_triang_getdSdA(double a)
{
	double alpha, alfMax, dPdA, r;

	// --- if a > area corresponding to sMax, then
	//     use slope between sFull & sMax
	alfMax = RECT_TRIANG_ALFMAX;
	alpha = a / m_aFull;
	if (alpha > alfMax)
		return (m_sFull - m_sMax) / ((1.0 - alfMax) * m_aFull);

	// --- use generic central difference method for very small a
	if (alpha <= 1.0e-30) return generic_getdSdA(a);

	// --- find deriv. of wetted perimeter
	if (a > m_aBot) dPdA = 2.0 / m_wMax;  // for upper rectangle
	else dPdA = m_rBot / sqrt(a * m_sBot);  // for triang. bottom

											// --- get hyd. radius & evaluate section factor derivative formula
	r = rect_triang_getRofA(a);
	return  (5. / 3. - (2. / 3.) * dPdA * r) * pow(r, 2. / 3.);
}

double xsectClass::rect_triang_getAofY(double y)
{
	if (y <= m_yBot) return y * y * m_sBot;         // below upper section
	else return m_aBot + (y - m_yBot) * m_wMax;  // above bottom section
}

double xsectClass::rect_triang_getWofY(double y)
{
	if (y <= m_yBot) return 2.0 * m_sBot * y;  // below upper section
	else return m_wMax;                               // above bottom section
}

double xsectClass::mod_basket_getYofA(double a)
{
	double alpha, y1;

	// --- water level below top of rectangular bottom
	if (a <= m_aFull - m_aBot) return a / m_wMax;

	// --- find unfilled top area / area of full circular top
	alpha = (m_aFull - a) / (PI * m_rBot * m_rBot);

	// --- find unfilled height
	if (alpha < 0.04) y1 = getYcircular(alpha);
	else                y1 = lookup(alpha, Y_Circ, N_Y_Circ);
	y1 = 2.0 * m_rBot * y1;

	// --- return difference between full height & unfilled height
	return m_yFull - y1;
}

double xsectClass::mod_basket_getRofA(double a)
{
	double y1, p, theta1;

	// --- water level is below top of rectangular bottom;
	//     return hyd. radius of rectangle
	if (a <= m_aFull - m_aBot)
		return a / (m_wMax + 2.0 * a / m_wMax);

	// --- find height of empty area
	y1 = m_yFull - mod_basket_getYofA(a);

	// --- find angle of circular arc corresponding to this height
	theta1 = 2.0 * acos(1.0 - y1 / m_rBot);

	// --- find perimeter of wetted portion of circular arc
	//     (angle of full circular opening was stored in sBot)
	p = (m_sBot - theta1) * m_rBot;

	// --- add on wetted perimeter of bottom rectangular area
	y1 = m_yFull - m_yBot;
	p = p + 2.0*y1 + m_wMax;

	// --- return area / wetted perimeter
	return a / p;
}

double xsectClass::mod_basket_getdSdA(double a)
{
	double r, dPdA;

	// --- if water level below top of rectangular bottom but not
	//     empty then use same code as for rectangular xsection
	if (a <= m_aFull - m_aBot && a / m_aFull > 1.0e-30)
	{
		r = a / (m_wMax + 2.0 * a / m_wMax);
		dPdA = 2.0 / m_wMax;
		return  (5. / 3. - (2. / 3.) * dPdA * r) * pow(r, 2. / 3.);
	}

	// --- otherwise use generic function
	else return generic_getdSdA(a);
}

double xsectClass::mod_basket_getAofY(double y)
{
	double a1, theta1, y1;

	// --- if water level is below top of rectangular bottom
	//     return depth * width
	if (y <= m_yFull - m_yBot) return y * m_wMax;

	// --- find empty top circular area
	y1 = m_yFull - y;
	theta1 = 2.0*acos(1.0 - y1 / m_rBot);
	a1 = 0.5 * m_rBot * m_rBot * (theta1 - sin(theta1));

	// --- return difference between full and empty areas
	return m_aFull - a1;
}

double xsectClass::mod_basket_getWofY(double y)
{
	double y1;

	// --- if water level below top of rectangular bottom then return width
	if (y <= 0.0) return 0.0;
	if (y <= m_yFull - m_yBot) return m_wMax;

	// --- find width of empty top circular section
	y1 = m_yFull - y;
	return 2.0 * sqrt(y1 * (2.0*m_rBot - y1));
}

double xsectClass::triang_getYofA(double a)
{
	return sqrt(a / m_sBot);
}

double xsectClass::triang_getRofA(double a)
{
	return a / (2. * triang_getYofA(a) * m_rBot);
}

double xsectClass::triang_getdSdA(double a)
{
	double r, dPdA;
	// --- use generic finite difference method for very small 'a'
	if (a / m_aFull <= 1.0e-30) return generic_getdSdA(a);

	// --- evaluate dSdA = [5/3 - (2/3)(dP/dA)R]R^(2/3)
	r = triang_getRofA(a);
	dPdA = m_rBot / sqrt(a * m_sBot);
	return  (5. / 3. - (2. / 3.) * dPdA * r) * pow(r, 2. / 3.);
}

double xsectClass::triang_getAofY(double y)
{
	return y * y * m_sBot;
}

double xsectClass::triang_getWofY(double y)
{
	return 2.0 * m_sBot * y;
}

double xsectClass::xsect_getAofY(double y)
//
//  Input:   xsect = ptr. to a cross section data structure
//           y = depth (ft)
//  Output:  returns area (ft2)
//  Purpose: computes xsection's area at a given depth.
//
{
	double yNorm = y / m_yFull;
	if (y <= 0.0) return 0.0;
	switch (m_type)
	{
	case FORCE_MAIN:
	case CIRCULAR:
		return m_aFull * lookup(yNorm, A_Circ, N_A_Circ);

	case FILLED_CIRCULAR:
		return filled_circ_getAofY(y);

	case EGGSHAPED:
		return m_aFull * lookup(yNorm, A_Egg, N_A_Egg);

	case HORSESHOE:
		return m_aFull * lookup(yNorm, A_Horseshoe, N_A_Horseshoe);

	case GOTHIC:
		return m_aFull * invLookup(yNorm, Y_Gothic, N_Y_Gothic);

	case CATENARY:
		return m_aFull * invLookup(yNorm, Y_Catenary, N_Y_Catenary);

	case SEMIELLIPTICAL:
		return m_aFull * invLookup(yNorm, Y_SemiEllip, N_Y_SemiEllip);

	case BASKETHANDLE:
		return m_aFull * lookup(yNorm, A_Baskethandle, N_A_Baskethandle);

	case SEMICIRCULAR:
		return m_aFull * invLookup(yNorm, Y_SemiCirc, N_Y_SemiCirc);

	case HORIZ_ELLIPSE:
		return m_aFull * lookup(yNorm, A_HorizEllipse, N_A_HorizEllipse);

	case VERT_ELLIPSE:
		return m_aFull * lookup(yNorm, A_VertEllipse, N_A_VertEllipse);

	case ARCH:
		return m_aFull * lookup(yNorm, A_Arch, N_A_Arch);

	case IRREGULAR:
		return m_aFull * lookup(yNorm,m_transectPtr->m_areaTbl, N_TRANSECT_TBL);

	case CUSTOM:
		return m_aFull * lookup(yNorm,m_shapePtr->m_areaTbl, N_SHAPE_TBL);

	case RECT_CLOSED:  return y * m_wMax;

	case RECT_TRIANG: return rect_triang_getAofY(y);

	case RECT_ROUND:  return rect_round_getAofY(y);

	case RECT_OPEN:   return y * m_wMax;

	case MOD_BASKET:  return mod_basket_getAofY(y);

	case TRAPEZOIDAL: return trapez_getAofY(y);

	case TRIANGULAR:  return triang_getAofY(y);

	case PARABOLIC:   return parab_getAofY(y);

	case POWERFUNC:   return powerfunc_getAofY(y);

	default:          return 0.0;
	}
}

double xsectClass::xsect_getRofY( double y)
//
//  Input:   xsect = ptr. to a cross section data structure
//           y = depth (ft)
//  Output:  returns hydraulic radius (ft)
//  Purpose: computes xsection's hydraulic radius at a given depth.
//
{
	double yNorm = y / m_yFull;
	switch (m_type)
	{
	case FORCE_MAIN:
	case CIRCULAR:
		return m_rFull * lookup(yNorm, R_Circ, N_R_Circ);

	case FILLED_CIRCULAR:
		if (m_yBot == 0.0)
			return m_rFull * lookup(yNorm, R_Circ, N_R_Circ);
		return filled_circ_getRofY(y);

	case EGGSHAPED:
		return m_rFull * lookup(yNorm, R_Egg, N_R_Egg);

	case HORSESHOE:
		return m_rFull * lookup(yNorm, R_Horseshoe, N_R_Horseshoe);

	case BASKETHANDLE:
		return m_rFull * lookup(yNorm, R_Baskethandle, N_R_Baskethandle);

	case HORIZ_ELLIPSE:
		return m_rFull * lookup(yNorm, R_HorizEllipse, N_R_HorizEllipse);

	case VERT_ELLIPSE:
		return m_rFull * lookup(yNorm, R_VertEllipse, N_R_VertEllipse);

	case ARCH:
		return m_rFull * lookup(yNorm, R_Arch, N_R_Arch);

	case IRREGULAR:
		return m_rFull * lookup(yNorm,m_transectPtr->m_hradTbl, N_TRANSECT_TBL);

	case CUSTOM:
		return m_rFull * lookup(yNorm,m_shapePtr->m_hradTbl, N_SHAPE_TBL);

	case RECT_TRIANG:  return rect_triang_getRofY(y);

	case RECT_ROUND:   return rect_round_getRofY(y);

	case TRAPEZOIDAL:  return trapez_getRofY(y);

	case TRIANGULAR:   return triang_getRofY(y);

	case PARABOLIC:    return parab_getRofY(y);

	case POWERFUNC:    return powerfunc_getRofY(y);

	default:           return xsect_getRofA(xsect_getAofY(y));
	}
}

double xsectClass::circ_getYofA(double a)
{
	double alpha = a / m_aFull;

	// --- use special function for small a/aFull
	if (alpha < 0.04)  return m_yFull * getYcircular(alpha);

	// --- otherwise use table
	else return m_yFull * lookup(alpha, Y_Circ, N_Y_Circ);
}

double xsectClass::getThetaOfPsi(double psi)
{
	int    k;
	double theta, theta1, ap, tt, tt23, t3, d;

	if (psi > 0.90)  theta = 4.17 + 1.12 * (psi - 0.90) / 0.176;
	else if (psi > 0.5)   theta = 3.14 + 1.03 * (psi - 0.5) / 0.4;
	else if (psi > 0.015) theta = 1.2 + 1.94 * (psi - 0.015) / 0.485;
	else                  theta = 0.12103 - 55.5075 * psi +
		15.62254 * sqrt(psi);
	theta1 = theta;
	ap = (2.0*PI) * psi;

	for (k = 1; k <= 40; k++)
	{
		theta = fabs(theta);
		tt = theta - sin(theta);
		tt23 = pow(tt, 2. / 3.);
		t3 = pow(theta, 1. / 3.);
		d = ap * theta / t3 - tt * tt23;
		d = d / (ap*(2. / 3.) / t3 - (5. / 3.)*tt23*(1.0 - cos(theta)));
		theta = theta - d;
		if (fabs(d) <= 0.0001) return theta;
	}
	return theta1;
}

double xsectClass::getAcircular(double psi)
{
	double theta;
	if (psi >= 1.0) return 1.0;
	if (psi <= 0.0) return 0.0;
	if (psi <= 1.0e-6)
	{
		theta = pow(124.4797*psi, 3. / 13.);
		return theta*theta*theta / 37.6911;
	}
	theta = getThetaOfPsi(psi);
	return (theta - sin(theta)) / (2.0 * PI);
}

double xsectClass::circ_getAofS(double s)
{
	double psi = s / m_sFull;
	if (psi == 0.0) return 0.0;
	if (psi >= 1.0) return m_aFull;

	// --- use special function for small s/sFull
	if (psi <= 0.015) return m_aFull * getAcircular(psi);

	// --- otherwise use table
	else return m_aFull * invLookup(psi, S_Circ, N_S_Circ);
}

double xsectClass::tabular_getdSdA(double a, double *table, int nItems)
//
//  Input:   xsect = ptr. to cross section data structure
//           a = area (ft2)
//           table = ptr. to table of section factor v. normalized area
//           nItems = number of equally spaced items in table
//  Output:  returns derivative of section factor w.r.t. area (ft^2/3)
//  Purpose: computes derivative of section factor w.r.t area
//           using geometry tables.
//
{
	int    i;
	double alpha = a / m_aFull;
	double delta = 1.0 / (nItems - 1);
	double dSdA;

	// --- find which segment of table contains alpha
	i = (int)(alpha / delta);
	if (i >= nItems - 1) i = nItems - 2;

	// --- compute slope from this interval of table
	dSdA = (table[i + 1] - table[i]) / delta;

	// --- convert slope to un-normalized value
	return dSdA * m_sFull / m_aFull;
}

double xsectClass::circ_getdSdA(double a)
{
	double alpha, theta, p, r, dPdA;

	// --- for near-zero area, use generic central difference formula
	alpha = a / m_aFull;
	if (alpha <= 1.0e-30) return 1.0e-30;  //generic_getdSdA(xsect, a);

										   // --- for small a/aFull use analytical derivative
	else if (alpha < 0.04)
	{
		theta = getThetaOfAlpha(alpha);
		p = theta * m_yFull / 2.0;
		r = a / p;
		dPdA = 4.0 / m_yFull / (1. - cos(theta));
		return  (5. / 3. - (2. / 3.) * dPdA * r) * pow(r, 2. / 3.);
	}

	// --- otherwise use generic tabular getdSdA
	else return tabular_getdSdA(a, S_Circ, N_S_Circ);
}

double xsectClass::filled_circ_getYofA(double a)
{
	double y;

	// --- temporarily remove filled portion of circle
	m_yFull += m_yBot;
	m_aFull += m_aBot;
	a += m_aBot;

	// --- find depth in unfilled circle
	y = circ_getYofA(a);

	// --- restore original values
	y -= m_yBot;
	m_yFull -= m_yBot;
	m_aFull -= m_aBot;
	return y;
}

double xsectClass::trapez_getYofA(double a)
{
	if (m_sBot == 0.0) return a / m_yBot;
	return (sqrt(m_yBot*m_yBot + 4.*m_sBot*a)
		- m_yBot) / (2. * m_sBot);
}

double xsectClass::trapez_getRofA(double a)
{
	return a / (m_yBot + trapez_getYofA(a) * m_rBot);
}

double xsectClass::trapez_getdSdA(double a)
{
	double r, dPdA;
	// --- use generic central difference method for very small a
	if (a / m_aFull <= 1.0e-30) return generic_getdSdA(a);

	// --- otherwise use analytical formula:
	//     dSdA = [5/3 - (2/3)(dP/dA)R]R^(2/3)
	r = trapez_getRofA(a);
	dPdA = m_rBot /
		sqrt(m_yBot * m_yBot + 4. * m_sBot * a);
	return  (5. / 3. - (2. / 3.) * dPdA * r) * pow(r, 2. / 3.);
}

double xsectClass::trapez_getWofY(double y)
{
	return m_yBot + 2.0 * y * m_sBot;
}

double xsectClass::xsect_getYofA(double a)
//
//  Input:   xsect = ptr. to a cross section data structure
//           a = area (ft2)
//  Output:  returns depth (ft)
//  Purpose: computes xsection's depth at a given area.
//
{
	double alpha = a / m_aFull;
	switch (m_type)
	{
	case FORCE_MAIN:
	case CIRCULAR: return circ_getYofA(a);

	case FILLED_CIRCULAR:
		return filled_circ_getYofA(a);

	case EGGSHAPED:
		return m_yFull * lookup(alpha, Y_Egg, N_Y_Egg);

	case HORSESHOE:
		return m_yFull * lookup(alpha, Y_Horseshoe, N_Y_Horseshoe);

	case GOTHIC:
		return m_yFull * lookup(alpha, Y_Gothic, N_Y_Gothic);

	case CATENARY:
		return m_yFull * lookup(alpha, Y_Catenary, N_Y_Catenary);

	case SEMIELLIPTICAL:
		return m_yFull * lookup(alpha, Y_SemiEllip, N_Y_SemiEllip);

	case BASKETHANDLE:
		return m_yFull * lookup(alpha, Y_BasketHandle, N_Y_BasketHandle);

	case SEMICIRCULAR:
		return m_yFull * lookup(alpha, Y_SemiCirc, N_Y_SemiCirc);

	case HORIZ_ELLIPSE:
		return m_yFull * invLookup(alpha, A_HorizEllipse, N_A_HorizEllipse);

	case VERT_ELLIPSE:
		return m_yFull * invLookup(alpha, A_VertEllipse, N_A_VertEllipse);

	case IRREGULAR:
		return m_yFull * invLookup(alpha,m_transectPtr->m_areaTbl, N_TRANSECT_TBL);

	case CUSTOM:
		return m_yFull * invLookup(alpha,m_shapePtr->m_areaTbl, N_SHAPE_TBL);

	case ARCH:
		return m_yFull * invLookup(alpha, A_Arch, N_A_Arch);

	case RECT_CLOSED: return a / m_wMax;

	case RECT_TRIANG: return rect_triang_getYofA(a);

	case RECT_ROUND:  return rect_round_getYofA(a);

	case RECT_OPEN:   return a / m_wMax;

	case MOD_BASKET:  return mod_basket_getYofA(a);

	case TRAPEZOIDAL: return trapez_getYofA(a);

	case TRIANGULAR:  return triang_getYofA(a);

	case PARABOLIC:   return parab_getYofA(a);

	case POWERFUNC:   return powerfunc_getYofA(a);

	default:          return 0.0;
	}
}

double xsectClass::xsect_getRofA(double a)
//
//  Input:   xsect = ptr. to a cross section data structure
//           a = area (ft2)
//  Output:  returns hydraulic radius (ft)
//  Purpose: computes xsection's hydraulic radius at a given area.
//
{
	double cathy;
	if (a <= 0.0) return 0.0;
	switch (m_type)
	{
	case HORIZ_ELLIPSE:
	case VERT_ELLIPSE:
	case ARCH:
	case IRREGULAR:
	case FILLED_CIRCULAR:
	case CUSTOM:
		return xsect_getRofY(xsect_getYofA(a));

	case RECT_CLOSED:  return rect_closed_getRofA(a);

	case RECT_OPEN:    return a / (m_wMax +
		(2. - m_sBot) * a / m_wMax);

	case RECT_TRIANG:  return rect_triang_getRofA(a);

	case RECT_ROUND:   return rect_round_getRofA(a);

	case MOD_BASKET:   return mod_basket_getRofA(a);

	case TRAPEZOIDAL:  return trapez_getRofA(a);

	case TRIANGULAR:   return triang_getRofA(a);

	case PARABOLIC:    return parab_getRofA(a);

	case POWERFUNC:    return powerfunc_getRofA(a);

	default:
		cathy = xsect_getSofA(a);
		if (cathy < TINY || a < TINY) return 0.0;
		return pow(cathy / a, 3. / 2.);
	}
}

double xsectClass::rect_closed_getSofA(double a)
{
	// --- if a > area corresponding to Smax then
	//     interpolate between sMax and Sfull
	double alfMax = RECT_ALFMAX;
	if (a / m_aFull > alfMax)
	{
		return m_sMax + (m_sFull - m_sMax) *
			(a / m_aFull - alfMax) / (1.0 - alfMax);
	}

	// --- otherwise use regular formula
	return a * pow(xsect_getRofA(a), 2. / 3.);
}

int xsectClass::xsect_setParams(int type, double p[], double ucf)
//
//  Input:   xsect = ptr. to a cross section data structure
//           type = xsection shape type
//           p[] = vector or xsection parameters
//           ucf = units correction factor
//  Output:  returns TRUE if successful, FALSE if not
//  Purpose: assigns parameters to a cross section's data structure.
//
{
    int    index;
    double aMax, theta;

    if ( type != DUMMY && p[0] <= 0.0 ) return FALSE;
    m_type  = type;
    switch ( type )
    {
    case DUMMY:
        m_yFull = TINY;
		m_wMax  = TINY;
		m_aFull = TINY;
		m_rFull = TINY;
		m_sFull = TINY;
		m_sMax  = TINY;
        break;

    case CIRCULAR:
		m_yFull = p[0]/ucf;
		m_wMax  = m_yFull;
		m_aFull = PI / 4.0 * m_yFull * m_yFull;
		m_rFull = 0.2500 * m_yFull;
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = 1.08 * m_sFull;
		m_ywMax = 0.5 * m_yFull;
        break;

    case FORCE_MAIN:
		m_yFull = p[0]/ucf;
		m_wMax  = m_yFull;
		m_aFull = PI / 4.0 * m_yFull * m_yFull;
		m_rFull = 0.2500 * m_yFull;
		m_sFull = m_aFull * pow(m_rFull, 0.63);
		m_sMax  = 1.06949 * m_sFull;
		m_ywMax = 0.5 * m_yFull;

        // --- save C-factor or roughness in rBot position
		m_rBot  = p[1];
        break;

    case FILLED_CIRCULAR:
        if ( p[1] >= p[0] ) return FALSE;

        // --- initially compute full values for unfilled pipe
		m_yFull = p[0]/ucf;
		m_wMax  = m_yFull;
		m_aFull = PI / 4.0 * m_yFull * m_yFull;
		m_rFull = 0.2500 * m_yFull;

        // --- find:
        //     yBot = depth of filled bottom
        //     aBot = area of filled bottom
        //     sBot = width of filled bottom
        //     rBot = wetted perimeter of filled bottom
		m_yBot  = p[1]/ucf;
		m_aBot  = circ_getAofY(m_yBot);
		m_sBot  = xsect_getWofY(m_yBot);
		m_rBot  = m_aBot / (m_rFull *
                       lookup(m_yBot/ m_yFull, R_Circ, N_R_Circ));

        // --- revise full values for filled bottom
		m_aFull -= m_aBot;
		m_rFull = m_aFull /
                       (PI*m_yFull - m_rBot + m_sBot);
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = 1.08 * m_sFull;
		m_yFull -= m_yBot;
		m_ywMax = 0.5 * m_yFull;
        break;

    case EGGSHAPED:
		m_yFull = p[0]/ucf;
		m_aFull = 0.5105 * m_yFull * m_yFull;
		m_rFull = 0.1931 * m_yFull;
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = 1.065 * m_sFull;
		m_wMax  = 2./3. * m_yFull;
		m_ywMax = 0.64 * m_yFull;
        break;

    case HORSESHOE:
		m_yFull = p[0]/ucf;
		m_aFull = 0.8293 * m_yFull * m_yFull;
		m_rFull = 0.2538 * m_yFull;
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = 1.077 * m_sFull;
		m_wMax  = 1.0 * m_yFull;
		m_ywMax = 0.5 * m_yFull;
        break;

    case GOTHIC:
		m_yFull = p[0]/ucf;
		m_aFull = 0.6554 * m_yFull * m_yFull;
		m_rFull = 0.2269 * m_yFull;
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = 1.065 * m_sFull;
		m_wMax  = 0.84 * m_yFull;
		m_ywMax = 0.45 * m_yFull;
        break;

    case CATENARY:
		m_yFull = p[0]/ucf;
		m_aFull = 0.70277 * m_yFull * m_yFull;
		m_rFull = 0.23172 * m_yFull;
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = 1.05 * m_sFull;
		m_wMax  = 0.9 * m_yFull;
		m_ywMax = 0.25 * m_yFull;
        break;

    case SEMIELLIPTICAL:
		m_yFull = p[0]/ucf;
		m_aFull = 0.785 * m_yFull * m_yFull;
		m_rFull = 0.242 * m_yFull;
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = 1.045 * m_sFull;
		m_wMax  = 1.0 * m_yFull;
		m_ywMax = 0.15 * m_yFull;
        break;

    case BASKETHANDLE:
		m_yFull = p[0]/ucf;
		m_aFull = 0.7862 * m_yFull * m_yFull;
		m_rFull = 0.2464 * m_yFull;
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = 1.06078 * m_sFull;
		m_wMax  = 0.944 * m_yFull;
		m_ywMax = 0.2 * m_yFull;
        break;

    case SEMICIRCULAR:
		m_yFull = p[0]/ucf;
		m_aFull = 1.2697 * m_yFull * m_yFull;
		m_rFull = 0.2946 * m_yFull;
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = 1.06637 * m_sFull;
		m_wMax  = 1.64 * m_yFull;
		m_ywMax = 0.15 * m_yFull;
        break;

    case RECT_CLOSED://p[0]是高，p[1]是宽，p[2]是淤积深度
        if ( p[1] <= 0.0 ) return FALSE;
		if (p[2] > 0 && p[2] < p[0])//输入了淤积深度的信息
		{
			m_yBot = p[2] / ucf;
			m_yFull = (p[0]-p[2]) / ucf;//矩形断面的高			
			m_wMax = p[1] / ucf;//矩形断面的宽
			m_aFull = m_yFull * m_wMax;////矩形断面的面积
			m_rFull = m_aFull / (2.0 * (m_yFull + m_wMax));//矩形断面满流时的水力半径
			m_sFull = m_aFull * pow(m_rFull, 2. / 3.);//矩形断面满流时A*R^(2/3)的值，A为满流时过流断面面积，R为满流时水力半径
			aMax = RECT_ALFMAX * m_aFull;//充满度0.97时的面积
			m_sMax = aMax * pow(rect_closed_getRofA(aMax), 2. / 3.);//充满度0.97时A*R^(2/3)的值，A、R分别为充满度0.97时过流断面面积与水力半径
			m_ywMax = m_yFull;//断面最宽处的宽度
			m_aBot = m_yBot*m_wMax;//淤积的面积
		}
		else//没有输入淤积深度，则按照原有的方式输入参数
		{
			m_yFull = p[0] / ucf;//矩形断面的高
			m_wMax = p[1] / ucf;//矩形断面的宽
			m_aFull = m_yFull * m_wMax;////矩形断面的面积
			m_rFull = m_aFull / (2.0 * (m_yFull + m_wMax));//矩形断面满流时的水力半径
			m_sFull = m_aFull * pow(m_rFull, 2. / 3.);//矩形断面满流时A*R^(2/3)的值，A为满流时过流断面面积，R为满流时水力半径
			aMax = RECT_ALFMAX * m_aFull;//充满度0.97时的面积
			m_sMax = aMax * pow(rect_closed_getRofA(aMax), 2. / 3.);//充满度0.97时A*R^(2/3)的值，A、R分别为充满度0.97时过流断面面积与水力半径
			m_ywMax = m_yFull;//断面最宽处的宽度
		}
        break;

    case RECT_OPEN:
        if ( p[1] <= 0.0 ) return FALSE;
		m_yFull = p[0]/ucf;
		m_wMax  = p[1]/ucf;
        if (p[2] < 0.0 || p[2] > 2.0) return FALSE;   //# sides to ignore
		m_sBot = p[2];
		m_aFull = m_yFull * m_wMax;
		m_rFull = m_aFull / ((2.0 - m_sBot) *
                       m_yFull + m_wMax);
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = m_sFull;
		m_ywMax = m_yFull;
        break;

    case RECT_TRIANG:
        if ( p[1] <= 0.0 || p[2] <= 0.0 ) return FALSE;
		m_yFull = p[0]/ucf;
		m_wMax  = p[1]/ucf;
		m_yBot  = p[2]/ucf;
		m_ywMax = m_yFull;

        // --- area of bottom triangle
		m_aBot  = m_yBot * m_wMax / 2.0;

        // --- slope of bottom side wall
		m_sBot  = m_wMax / m_yBot / 2.0;

        // --- length of side wall per unit of depth
		m_rBot  = sqrt( 1. + m_sBot * m_sBot );

		m_aFull = m_wMax * (m_yFull - m_yBot / 2.0);
		m_rFull = m_aFull / (2.0 * m_yBot * m_rBot + 2.0 *
                        (m_yFull - m_yBot) + m_wMax);
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
        aMax = RECT_TRIANG_ALFMAX * m_aFull;
		m_sMax  = aMax * pow(rect_triang_getRofA(aMax), 2./3.);
        break;

    case RECT_ROUND:
        if ( p[1] <= 0.0 ) return FALSE;
        if ( p[2] < p[1]/2.0 ) p[2] = p[1]/2.0;
		m_yFull = p[0]/ucf;
		m_wMax  = p[1]/ucf;
		m_rBot  = p[2]/ucf;

        // --- angle of circular arc
        theta = 2.0 * asin(m_wMax / 2.0 / m_rBot);

        // --- area of circular bottom
		m_aBot  = m_rBot * m_rBot /
                       2.0 * (theta - sin(theta));

        // --- section factor for circular bottom
		m_sBot  = PI * m_rBot * m_rBot *
                       pow(m_rBot/2.0, 2./3.);

        // --- depth of circular bottom
		m_yBot  = m_rBot * (1.0 - cos(theta/2.0));
		m_ywMax = m_yFull;

		m_aFull = m_wMax * (m_yFull - m_yBot) + m_aBot;
		m_rFull = m_aFull / (m_rBot * theta + 2.0 *
                        (m_yFull - m_yBot) + m_wMax);
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
        aMax = RECT_ROUND_ALFMAX * m_aFull;
		m_sMax = aMax * pow(rect_round_getRofA(aMax), 2./3.);
        break;

    case MOD_BASKET:
        if ( p[1] <= 0.0 ) return FALSE;
        if ( p[2] < p[1]/2.0 ) p[2] = p[1]/2.0;
		m_yFull = p[0]/ucf;
		m_wMax  = p[1]/ucf;

        // --- radius of circular arc
		m_rBot = p[2]/ucf;

        // --- angle of circular arc
        theta = 2.0 * asin(m_wMax / 2.0 / m_rBot);
		m_sBot = theta;

        // --- height of circular arc
		m_yBot = m_rBot * (1.0 - cos(theta/2.0));
		m_ywMax = m_yFull - m_yBot;

        // --- area of circular arc
		m_aBot = m_rBot * m_rBot /
                      2.0 * (theta - sin(theta));

        // --- full area
		m_aFull = (m_yFull - m_yBot) * m_wMax +
			m_aBot;

        // --- full hydraulic radius & section factor
		m_rFull = m_aFull / (m_rBot * theta + 2.0 *
                        (m_yFull - m_yBot) + m_wMax);
		m_sFull = m_aFull * pow(m_rFull, 2./3.);

        // --- area corresponding to max. section factor
		m_sMax = xsect_getSofA(Amax[MOD_BASKET]* m_aFull);
        break;

    case TRAPEZOIDAL:
        if ( p[1] < 0.0 || p[2] < 0.0 || p[3] < 0.0 ) return FALSE;
		m_yFull = p[0]/ucf;
		m_ywMax = m_yFull;

        // --- bottom width
		m_yBot = p[1]/ucf;
		//added by jinxi
		m_rightSlope = p[2];//geom3
		m_leftSlope = p[3];//geom4
								//end of added by jinxi
        // --- avg. slope of side walls
		m_sBot  = ( p[2] + p[3] )/2.0;
        if (m_yBot == 0.0 && m_sBot == 0.0 ) return FALSE;

        // --- length of side walls per unit of depth
		m_rBot  = sqrt( 1.0 + p[2]*p[2] ) + sqrt( 1.0 + p[3]*p[3] );

        // --- top width
		m_wMax = m_yBot + m_yFull * (p[2] + p[3]);

		m_aFull = (m_yBot + m_sBot * m_yFull ) * m_yFull;
		m_rFull = m_aFull / (m_yBot + m_yFull * m_rBot);
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = m_sFull;
        break;

    case TRIANGULAR:
        if ( p[1] <= 0.0 ) return FALSE;
		m_yFull = p[0]/ucf;
		m_wMax  = p[1]/ucf;
		m_ywMax = m_yFull;

        // --- slope of side walls
		m_sBot  = m_wMax / m_yFull / 2.;

        // --- length of side wall per unit of depth
		m_rBot  = sqrt( 1. + m_sBot * m_sBot );

		m_aFull = m_yFull * m_yFull * m_sBot;
		m_rFull = m_aFull / (2.0 * m_yFull * m_rBot);
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = m_sFull;
        break;

    case PARABOLIC:
        if ( p[1] <= 0.0 ) return FALSE;
		m_yFull = p[0]/ucf;
		m_wMax  = p[1]/ucf;
		m_ywMax = m_yFull;

        // --- rBot :: 1/c^.5, where y = c*x^2 is eqn. of parabolic shape
		m_rBot  = m_wMax / 2.0 / sqrt(m_yFull);

		m_aFull = (2./3.) * m_yFull * m_wMax;
		m_rFull = xsect_getRofY(m_yFull);
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = m_sFull;
        break;

    case POWERFUNC:
        if ( p[1] <= 0.0 || p[2] <= 0.0 ) return FALSE;
		m_yFull = p[0]/ucf;
		m_wMax  = p[1]/ucf;
		m_ywMax = m_yFull;
		m_sBot  = 1.0 / p[2];
		m_rBot  = m_wMax / (m_sBot + 1) /
                       pow(m_yFull, m_sBot);
		m_aFull = m_yFull * m_wMax / (m_sBot+1);
		m_rFull = xsect_getRofY(m_yFull);
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = m_sFull;
        break;

    case HORIZ_ELLIPSE:
        if ( p[1] == 0.0 ) p[2] = p[0];
        if ( p[2] > 0.0 )                        // std. ellipse pipe
        {
            index = (int)floor(p[2]) - 1;        // size code
            if ( index < 0 ||
                 index >= NumCodesEllipse ) return FALSE;
			m_yFull = MinorAxis_Ellipse[index]/12.;
			m_wMax  = MajorAxis_Ellipse[index]/12.;
			m_aFull = Afull_Ellipse[index];
			m_rFull = Rfull_Ellipse[index];
        }
        else
        {
            // --- length of minor axis
			m_yFull = p[0]/ucf;

            // --- length of major axis
            if ( p[1] < 0.0 ) return FALSE;
			m_wMax = p[1]/ucf;
			m_aFull = 1.2692 * m_yFull * m_yFull;
			m_rFull = 0.3061 * m_yFull;
        }
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = m_sFull;
		m_ywMax = 0.48 * m_yFull;
        break;

    case VERT_ELLIPSE:
        if ( p[1] == 0.0 ) p[2] = p[0];
        if ( p[2] > 0.0 )                        // std. ellipse pipe
        {
            index = (int)floor(p[2]) - 1;        // size code
            if ( index < 0 ||
                 index >= NumCodesEllipse ) return FALSE;
			m_yFull = MajorAxis_Ellipse[index]/12.;
			m_wMax  = MinorAxis_Ellipse[index]/12.;
			m_aFull = Afull_Ellipse[index];
			m_rFull = Rfull_Ellipse[index];
        }
        else
        {
            // --- length of major axis
            if ( p[1] < 0.0 ) return FALSE;

            // --- length of minor axis
			m_yFull = p[0]/ucf;
			m_wMax = p[1]/ucf;
			m_aFull = 1.2692 * m_wMax * m_wMax;
			m_rFull = 0.3061 * m_wMax;
        }
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = m_sFull;
		m_ywMax = 0.48 * m_yFull;
        break;

    case ARCH:
        if ( p[1] == 0.0 ) p[2] = p[0];
        if ( p[2] > 0.0 )                        // std. arch pipe
        {
            index = (int)floor(p[2]) - 1;        // size code
            if ( index < 0 ||
                 index >= NumCodesArch ) return FALSE;
			m_yFull = Yfull_Arch[index]/12.;     // Yfull units are inches
			m_wMax  = Wmax_Arch[index]/12.;      // Wmax units are inches
			m_aFull = Afull_Arch[index];
			m_rFull = Rfull_Arch[index];
        }
        else                                     // non-std. arch pipe
        {
            if ( p[1] < 0.0 ) return FALSE;
			m_yFull = p[0]/ucf;
			m_wMax  = p[1]/ucf;
			m_aFull = 0.7879 * m_yFull * m_wMax;
			m_rFull = 0.2991 * m_yFull;
        }
		m_sFull = m_aFull * pow(m_rFull, 2./3.);
		m_sMax  = m_sFull;
		m_ywMax = 0.28 * m_yFull;
        break;
    }
    return TRUE;
}

void xsectClass::xsect_setIrregXsectParams()
//
//  Input:   xsect = ptr. to a cross section data structure
//  Output:  none
//  Purpose: assigns transect parameters to an irregular shaped cross section.
//
{
    int index = m_transect;
    int     i, iMax;
    double  wMax;
    double* wTbl = m_transectPtr->m_widthTbl;

    m_yFull = m_transectPtr->m_yFull;
	m_wMax  = m_transectPtr->m_wMax;
	m_aFull = m_transectPtr->m_aFull;
	m_rFull = m_transectPtr->m_rFull;
	m_sFull = m_aFull * pow(m_rFull, 2./3.);
	m_sMax = m_transectPtr->m_sMax;
	m_aBot = m_transectPtr->m_aMax;

    // Search transect's width table up to point where width decreases
    iMax = 0;
    wMax = wTbl[0];
    for (i = 1; i < N_TRANSECT_TBL; i++)
    {
	if ( wTbl[i] < wMax ) break;
	wMax = wTbl[i];
	iMax = i;
    }

    // Determine height at lowest widest point
    m_ywMax = m_yFull * (double)iMax / (double)(N_TRANSECT_TBL-1);
}

void xsectClass::xsect_setCustomXsectParams()
//
//  Input:   xsect = ptr. to a cross section data structure
//  Output:  none
//  Purpose: assigns parameters to a custom-shaped cross section.
//
{
    double  yFull = m_yFull;
    int     i, iMax;
    double  wMax;
    double* wTbl =  m_shapePtr->m_widthTbl;

    m_wMax  = m_shapePtr->m_wMax * yFull;
    m_aFull = m_shapePtr->m_aFull * yFull * yFull;
	m_rFull = m_shapePtr->m_rFull * yFull;
	m_sFull = m_aFull * pow(m_rFull, 2./3.);
	m_sMax  = m_shapePtr->m_sMax * yFull * yFull * pow(yFull, 2./3.);
	m_aBot  = m_shapePtr->m_aMax * yFull * yFull;

    // Search shape's width table up to point where width decreases
    iMax = 0;
    wMax = wTbl[0];
    for (i = 1; i < N_SHAPE_TBL; i++)
    {
	if ( wTbl[i] < wMax ) break;
	wMax = wTbl[i];
	iMax = i;
    }

    // Determine height at lowest widest point
    m_ywMax = yFull * (double)iMax / (double)(N_SHAPE_TBL-1);
}

double xsectClass::xsect_getAmax()
//
//  Input:   xsect = ptr. to a cross section data structure
//  Output:  returns area (ft2)
//  Purpose: finds xsection area at maximum flow depth.
//
{
    if ( m_type == IRREGULAR ) return m_aBot;
    else if (m_type == CUSTOM ) return m_aBot;
    else return Amax[m_type] * m_aFull;
}

double xsectClass::rect_round_getdSdA(double a)
{
	double alfMax, r, dPdA;

	// --- if a > area corresponding to sMax, then
	//     use slope between sFull & sMax
	alfMax = RECT_ROUND_ALFMAX;
	if (a / m_aFull > alfMax)
	{
		return (m_sFull - m_sMax) /
			((1.0 - alfMax) * m_aFull);
	}

	// --- if above circular invert, use analytical function for dS/dA
	else if (a > m_aBot)
	{
		r = rect_round_getRofA(a);
		dPdA = 2.0 / m_wMax;       // d(wet perim)/dA for rect.
		return  (5. / 3. - (2. / 3.) * dPdA * r) * pow(r, 2. / 3.);
	}

	// --- otherwise use generic finite difference function
	else return generic_getdSdA(a);
}

double xsectClass::rect_round_getWofY(double y)
{
	// --- return width if depth above circular bottom section
	if (y > m_yBot) return m_wMax;

	// --- find width of circular section
	return 2.0 * sqrt(y * (2.0*m_rBot - y));
}

double xsectClass::xsect_getWofY(double y)
//
//  Input:   xsect = ptr. to a cross section data structure
//           y = depth ft)
//  Output:  returns top width (ft)
//  Purpose: computes xsection's top width at a given depth.
//
{
    double yNorm = y / m_yFull;
    switch ( m_type )
    {
      case FORCE_MAIN:
      case CIRCULAR:
        return m_wMax * lookup(yNorm, W_Circ, N_W_Circ);

      case FILLED_CIRCULAR:
        yNorm = (y + m_yBot) / (m_yFull + m_yBot);
        return m_wMax * lookup(yNorm, W_Circ, N_W_Circ);

      case EGGSHAPED:
        return m_wMax * lookup(yNorm, W_Egg, N_W_Egg);

      case HORSESHOE:
        return m_wMax * lookup(yNorm, W_Horseshoe, N_W_Horseshoe);

      case GOTHIC:
        return m_wMax * lookup(yNorm, W_Gothic, N_W_Gothic);

      case CATENARY:
        return m_wMax * lookup(yNorm, W_Catenary, N_W_Catenary);

      case SEMIELLIPTICAL:
        return m_wMax * lookup(yNorm, W_SemiEllip, N_W_SemiEllip);

      case BASKETHANDLE:
        return m_wMax * lookup(yNorm, W_BasketHandle, N_W_BasketHandle);

      case SEMICIRCULAR:
        return m_wMax * lookup(yNorm, W_SemiCirc, N_W_SemiCirc);

      case HORIZ_ELLIPSE:
        return m_wMax * lookup(yNorm, W_HorizEllipse, N_W_HorizEllipse);

      case VERT_ELLIPSE:
        return m_wMax * lookup(yNorm, W_VertEllipse, N_W_VertEllipse);

      case ARCH:
        return m_wMax * lookup(yNorm, W_Arch, N_W_Arch);

      case IRREGULAR:
        return m_wMax * lookup(yNorm,m_transectPtr->m_widthTbl, N_TRANSECT_TBL);

      case CUSTOM:
        return m_wMax * lookup(yNorm,m_shapePtr->m_widthTbl, N_SHAPE_TBL);

      case RECT_CLOSED: 
          if (yNorm == 1.0) return 0.0;                                        //(5.1.013)
          return m_wMax;

      case RECT_TRIANG: return rect_triang_getWofY(y);

      case RECT_ROUND:  return rect_round_getWofY(y);

      case RECT_OPEN:   return m_wMax;

      case MOD_BASKET:  return mod_basket_getWofY(y);

      case TRAPEZOIDAL: return trapez_getWofY(y);

      case TRIANGULAR:  return triang_getWofY(y);

      case PARABOLIC:   return parab_getWofY(y);

      case POWERFUNC:   return powerfunc_getWofY(y);

      default:          return 0.0;
    }
}

double xsectClass::rect_closed_getdSdA(double a)
{
	double alpha, alfMax, r;

	// --- if above level corresponding to sMax, then
	//     use slope between sFull & sMax
	alfMax = RECT_ALFMAX;
	alpha = a / m_aFull;
	if (alpha > alfMax)
	{
		return (m_sFull - m_sMax) /
			((1.0 - alfMax) * m_aFull);
	}

	// --- for small a/aFull use generic central difference formula
	if (alpha <= 1.0e-30) return generic_getdSdA(a);

	// --- otherwise evaluate dSdA = [5/3 - (2/3)(dP/dA)R]R^(2/3)
	//     (where P = wetted perimeter & dPdA = 2/width)
	r = xsect_getRofA(a);
	return  (5. / 3. - (2. / 3.) * (2.0 / m_wMax) * r) * pow(r, 2. / 3.);
}

double xsectClass::rect_open_getdSdA(double a)
{
	double r, dPdA;

	// --- for small a/aFull use generic central difference formula
	if (a / m_aFull <= 1.0e-30) return generic_getdSdA(a);

	// --- otherwise evaluate dSdA = [5/3 - (2/3)(dP/dA)R]R^(2/3)
	//     (where P = wetted perimeter)
	r = xsect_getRofA(a);
	dPdA = (2.0 - m_sBot) / m_wMax; // since P = geom2 + 2a/geom2
	return  (5. / 3. - (2. / 3.) * dPdA * r) * pow(r, 2. / 3.);
}

double xsectClass::xsect_getdSdA(double a)
//
//  Input:   xsect = ptr. to a cross section data structure
//           a = area (ft2)
//  Output:  returns derivative of section factor w.r.t. area (ft^2/3)
//  Purpose: computes xsection's derivative of its section factor with
//           respect to area at a given area.
//
{
	switch (m_type)
	{
	case FORCE_MAIN:
	case CIRCULAR:
		return circ_getdSdA(a);

	case EGGSHAPED:
		return tabular_getdSdA(a, S_Egg, N_S_Egg);

	case HORSESHOE:
		return tabular_getdSdA(a, S_Horseshoe, N_S_Horseshoe);

	case GOTHIC:
		return tabular_getdSdA(a, S_Gothic, N_S_Gothic);

	case CATENARY:
		return tabular_getdSdA(a, S_Catenary, N_S_Catenary);

	case SEMIELLIPTICAL:
		return  tabular_getdSdA(a, S_SemiEllip, N_S_SemiEllip);

	case BASKETHANDLE:
		return  tabular_getdSdA(a, S_BasketHandle, N_S_BasketHandle);

	case SEMICIRCULAR:
		return  tabular_getdSdA(a, S_SemiCirc, N_S_SemiCirc);

	case RECT_CLOSED:
		return rect_closed_getdSdA(a);

	case RECT_OPEN:
		return rect_open_getdSdA(a);

	case RECT_TRIANG:
		return rect_triang_getdSdA(a);

	case RECT_ROUND:
		return rect_round_getdSdA(a);

	case MOD_BASKET:
		return mod_basket_getdSdA(a);

	case TRAPEZOIDAL:
		return trapez_getdSdA(a);

	case TRIANGULAR:
		return triang_getdSdA(a);

	default: return generic_getdSdA(a);
	}
}

void xsectClass::evalSofA(double a, double* f, double* df, void* p)
//
//  Input:   a = area
//  Output:  f = root finding function
//           df = derivative of root finding function
//  Purpose: function used in conjunction with getAofS() that evaluates
//           f = S(a) - s and df = dS(a)/dA.
//
{
	xsectStarClass* xsectStar;
	double s;

	xsectStar = (xsectStarClass *)p;
	s = xsectStar->m_xsect->xsect_getSofA(a);
	*f = s - xsectStar->m_s;
	*df = xsectStar->m_xsect->xsect_getdSdA(a);
}


double xsectClass::generic_getAofS(double s)
//
//  Input:   xsect = ptr. to a cross section data structure
//           s = section factor (ft^8/3)，s是Q/(n*√i)的计算结果
//  Output:  returns area (ft2)
//  Purpose: finds area given section factor by
//           solving S = A*(A/P(A))^(2/3) using Newton-Raphson iterations.这里的P(A)是湿周的意思
//
{
	double a, a1, a2, tol;
	xsectStarClass xsectStar;
	if (s <= 0.0) return 0.0;

	// --- if S is between sMax and sFull then
	//     bracket A between aFull and aMax
	if ((s <= m_sMax && s >= m_sFull)
		&& m_sMax != m_sFull)
	{
		a1 = m_aFull;          // do this because sFull < sMax
		a2 = xsect_getAmax();
	}

	// --- otherwise bracket A between 0 and aMax
	else
	{
		a1 = 0.0;
		a2 = xsect_getAmax();
	}

	// --- place S & xsect in xsectStar for access by evalSofA function
	xsectStar.m_xsect = this;
	xsectStar.m_s = s;

	// --- compute starting guess for A
	a = 0.5 * (a1 + a2);

	// use the Newton-Raphson root finder function to find A
	tol = 0.0001 * m_aFull;
	findroot_Newton(a1, a2, &a, tol, evalSofA, &xsectStar);
	return a;
}

double xsectClass::getQcritical(double yc, void* p)
//
//  Input:   yc = critical depth (ft)
//           p = pointer to a TXsectStar object
//  Output:  returns flow difference value (cfs)
//  Purpose: finds difference between critical flow at depth yc and
//           some target value.
//
{
	double a, w, qc;
	xsectStarClass* xsectStar;

	xsectStar = (xsectStarClass *)p;
	a = xsectStar->m_xsect->xsect_getAofY(yc);
	w = xsectStar->m_xsect->xsect_getWofY(yc);
	qc = -xsectStar->m_qc;
	if (w > 0.0)  qc = a * sqrt(GRAVITY * a / w) - xsectStar->m_qc;
	return qc;
}

double xsectClass::getYcritEnum(double q, double y0)
//
//  Input:   xsect = ptr. to cross section data structure
//           q = critical flow rate (cfs)
//           y0 = estimate of critical depth (ft)
//  Output:  returns true critical depth (ft)
//  Purpose: solves a * sqrt(a(y)*g / w(y)) - q for y using interval
//           enumeration with starting guess of y0.
//
{
	double     q0, dy, qc, yc;
	int        i1, i;
	xsectStarClass xsectStar;

	// --- divide cross section depth into 25 increments and
	//     locate increment corresponding to initial guess y0
	dy = m_yFull / 25.;
	i1 = (int)(y0 / dy);

	// --- evaluate critical flow at this increment
	xsectStar.m_xsect = this;
	xsectStar.m_qc = 0.0;
	q0 = getQcritical(i1*dy, &xsectStar);

	// --- initial flow lies below target flow
	if (q0 < q)
	{
		// --- search each successive higher depth increment
		yc = m_yFull;
		for (i = i1 + 1; i <= 25; i++)
		{
			// --- if critical flow at current depth is above target
			//     then use linear interpolation to compute critical depth
			qc = getQcritical(i*dy, &xsectStar);
			if (qc >= q)
			{
				yc = ((q - q0) / (qc - q0) + (double)(i - 1)) * dy;
				break;
			}
			q0 = qc;
		}
	}

	// --- initial flow lies above target flow
	else
	{
		// --- search each successively lower depth increment
		yc = 0.0;
		for (i = i1 - 1; i >= 0; i--)
		{
			// --- if critical flow at current depth is below target
			//     then use linear interpolation to compute critical depth
			qc = getQcritical(i*dy, &xsectStar);
			if (qc < q)
			{
				yc = ((q - qc) / (q0 - qc) + (double)i) * dy;
				break;
			}
			q0 = qc;
		}
	}
	return yc;
}

double xsectClass::getYcritRidder(double q, double y0)
//
//  Input:   xsect = ptr. to cross section data structure
//           q = critical flow rate (cfs)
//           y0 = estimate of critical depth (ft)
//  Output:  returns true critical depth (ft)
//  Purpose: solves a * sqrt(a(y)*g / w(y)) - q for y using Ridder's
//           root finding method with starting guess of y0.
//
{
	double  y1 = 0.0;
	double  y2 = 0.99 * m_yFull;
	double  yc;
	double q0, q1, q2;
	xsectStarClass xsectStar;

	// --- store reference to cross section in global pointer
	xsectStar.m_xsect = this;
	xsectStar.m_qc = 0.0;

	// --- check if critical flow at (nearly) full depth < target flow
	q2 = getQcritical(y2, &xsectStar);
	if (q2 < q) return m_yFull;

	// --- evaluate critical flow at initial depth guess y0
	//     and at 1/2 of full depth
	q0 = getQcritical(y0, &xsectStar);
	q1 = getQcritical(0.5*m_yFull, &xsectStar);

	// --- adjust search interval on depth so it contains flow q
	if (q0 > q)
	{
		y2 = y0;
		if (q1 < q) y1 = 0.5*m_yFull;
	}
	else
	{
		y1 = y0;
		if (q1 > q) y2 = 0.5*m_yFull;
	}

	// --- save value of target critical flow in global variable
	xsectStar.m_qc = q;

	// --- call Ridder root finding procedure with error tolerance
	//     of 0.001 ft. to find critical depth yc
	yc = findroot_Ridder(y1, y2, 0.001, getQcritical, &xsectStar);
	return yc;
}

double xsectClass::xsect_getAofS(double s)
//
//  Input:   xsect = ptr. to a cross section data structure
//           s = section factor (ft^(8/3))
//  Output:  returns area (ft2)
//  Purpose: computes xsection's area at a given section factor.
//
{
    double psi = s / m_sFull;
    if ( s <= 0.0 ) return 0.0;
    if ( s > m_sMax ) s = m_sMax;
    switch ( m_type )
    {
      case DUMMY:     return 0.0;

      case FORCE_MAIN:
      case CIRCULAR:  return circ_getAofS(s);

      case EGGSHAPED:
        return m_aFull * invLookup(psi, S_Egg, N_S_Egg);

      case HORSESHOE:
        return m_aFull * invLookup(psi, S_Horseshoe, N_S_Horseshoe);

      case GOTHIC:
        return m_aFull * invLookup(psi, S_Gothic, N_S_Gothic);

      case CATENARY:
        return m_aFull * invLookup(psi, S_Catenary, N_S_Catenary);

      case SEMIELLIPTICAL:
        return m_aFull * invLookup(psi, S_SemiEllip, N_S_SemiEllip);

      case BASKETHANDLE:
        return m_aFull * invLookup(psi, S_BasketHandle, N_S_BasketHandle);

      case SEMICIRCULAR:
        return m_aFull * invLookup(psi, S_SemiCirc, N_S_SemiCirc);

      default: return generic_getAofS(s);
    }
}



//=============================================================================

double xsectClass::xsect_getYcrit(double q)
//
//  Input:   xsect = ptr. to a cross section data structure
//           q = flow rate (cfs)
//  Output:  returns critical depth (ft)
//  Purpose: computes critical depth at a specific flow rate.
//
{
    double q2g = SQR(q) / GRAVITY;
    double y, r;

    if ( q2g == 0.0 ) return 0.0;
    switch ( m_type )
    {
      case DUMMY:
        return 0.0;

      case RECT_OPEN:
      case RECT_CLOSED:
        // --- analytical expression for yCritical is
        //     y = (q2g / w^2)^(1/3) where w = width
        y = pow(q2g / SQR(m_wMax), 1./3.);
        break;

      case TRIANGULAR:
        // --- analytical expression for yCritical is
        //     y = (2 * q2g / s^2)^(1/5) where s = side slope
        y = pow(2.0 * q2g / SQR(m_sBot), 1./5.);
        break;

      case PARABOLIC:
        // --- analytical expression for yCritical is
        //     y = (27/32 * q2g * c)^(1/4) where y = c*x^2
        //     is eqn. for parabola and 1/sqrt(c) = rBot
        y = pow(27./32. * q2g / SQR(m_rBot), 1./4.);
        break;

      case POWERFUNC:
        y = 1. / (2.0 * m_sBot + 3.0);
        y = pow( q2g * (m_sBot + 1.0) / SQR(m_rBot), y);
        break;

      default:
        // --- first estimate yCritical for an equivalent circular conduit
        //     using 1.01 * (q2g / yFull)^(1/4)
        y = 1.01 * pow(q2g / m_yFull, 1./4.);
        if (y >= m_yFull) y = 0.97 * m_yFull;

        // --- then find ratio of conduit area to equiv. circular area
        r = m_aFull / (PI / 4.0 * SQR(m_yFull));

        // --- use interval enumeration method to find yCritical if
        //     area ratio not too far from 1.0
        if ( r >= 0.5 && r <= 2.0 )
            y = getYcritEnum(q, y);

        // --- otherwise use Ridder's root finding method
        else y = getYcritRidder(q, y);
    }

    // --- do not allow yCritical to be > yFull
    return MIN(y, m_yFull);
}

//=============================================================================



//=============================================================================













double xsectClass::circ_getAofY(double y)
{
    double yNorm;
    yNorm = y / m_yFull;
    return m_aFull * lookup(yNorm, A_Circ, N_A_Circ);
}


#pragma region of jinxi
double xsectClass::linkScour_get_k(double x, int nItems)
{
	int itemIndex = locate(x, X_k_m, nItems - 1);//x的值在X_k_m[itemIndex]与X_k_m[itemIndex+1]之间
	if (itemIndex == nItems - 1)//传入的x值大于等于table中最大的一个值
	{
		return R_k[nItems - 1];
	}
	else
	{
		double x0 = X_k_m[itemIndex];
		double x1 = X_k_m[itemIndex + 1];
		double y0 = R_k[itemIndex];
		double y1 = R_k[itemIndex + 1];
		return y0 + (x - x0)*(y1 - y0) / (x1 - x0);
	}
}

double xsectClass::linkScour_get_m(double x, int nItems)
{
	int itemIndex = locate(x, X_k_m, nItems - 1);//x的值在X_k_m[itemIndex]与X_k_m[itemIndex+1]之间
	if (itemIndex == nItems - 1)//传入的x值大于等于table中最大的一个值
	{
		return R_m[nItems - 1];
	}
	else
	{
		double x0 = X_k_m[itemIndex];
		double x1 = X_k_m[itemIndex + 1];
		double y0 = R_m[itemIndex];
		double y1 = R_m[itemIndex + 1];
		return y0 + (x - x0)*(y1 - y0) / (x1 - x0);
	}
}

double xsectClass::linkScour_get_V(double x, int nItems)
{
	int itemIndex = locate(x, sedimentP_D, nItems - 1);//x的值在sedimentP_D[itemIndex]与X_k_m[itemIndex+1]之间,找到当前粒径的位置
	if (itemIndex == nItems - 1)//传入的x值大于等于table中最大的一个值
	{
		return sedimentP_V[nItems - 1];
	}
	else
	{
		double x0 = sedimentP_D[itemIndex];
		double x1 = sedimentP_D[itemIndex + 1];
		double y0 = sedimentP_V[itemIndex];
		double y1 = sedimentP_V[itemIndex + 1];
		return y0 + (x - x0)*(y1 - y0) / (x1 - x0);
	}
}

//冲淤计算过程中管底淤积深度会发生变化，通过该函数更新管道断面面积
void    projectClass::xsect_updateCrossSection(int linkIndex)
{
	double aMax;
	//根据淤积深度更新断面参数
	xsectClass *xsect = &(Link[linkIndex].m_xsect);
	switch (xsect->m_type)
	{
	case FILLED_CIRCULAR:
		//未淤积情况下的各项参数要按照未淤积的情况赋值yFull、wMax、aFull、rFull，这些参数会根据淤积情况修正，所以这些参数保存的是淤积状态下的全深、全宽、总面积、
		//满流水力半径，在更新淤积深度时这些参数要恢复到未淤积状态时的值
		xsect->m_yFull = xsect->m_wMax;//wMax在根据淤积深度更新过流断面参数时值没有变，而且对于圆管该值就是直径
		xsect->m_aFull = PI / 4.0 * xsect->m_yFull * xsect->m_yFull;
		xsect->m_rFull = 0.2500 * xsect->m_yFull;
		// --- 更新下列参数:
		//     yBot = depth of filled bottom
		//     aBot = area of filled bottom
		//     sBot = width of filled bottom
		//     rBot = wetted perimeter of filled bottom
		xsect->m_yBot = Link[linkIndex].m_sedimentH;
		xsect->m_aBot = xsect->circ_getAofY(xsect->m_yBot);
		xsect->m_sBot = xsect->xsect_getWofY(xsect->m_yBot);
		xsect->m_rBot = xsect->m_aBot / (xsect->m_rFull *	xsectClass::lookup(xsect->m_yBot / xsect->m_yFull, R_Circ, N_R_Circ));

		// --- revise full values for filled bottom
		xsect->m_aFull -= xsect->m_aBot;
		xsect->m_rFull = xsect->m_aFull / (PI*xsect->m_yFull - xsect->m_rBot + xsect->m_sBot);
		xsect->m_sFull = xsect->m_aFull * pow(xsect->m_rFull, 2. / 3.);
		xsect->m_sMax = 1.08 * xsect->m_sFull;
		xsect->m_yFull -= xsect->m_yBot;
		xsect->m_ywMax = 0.5 * xsect->m_yFull;
		//根据淤积深度更新管段上下游端点处的offset
		Link[linkIndex].m_offset1 += Link[linkIndex].m_sedimentDeltaH;
		Link[linkIndex].m_offset2 += Link[linkIndex].m_sedimentDeltaH;
		break;
	case RECT_CLOSED:
		xsect->m_yFull = xsect->m_yBot + xsect->m_yFull - Link[linkIndex].m_sedimentH;//矩形断面的原始高，不考虑淤积
		xsect->m_yBot = Link[linkIndex].m_sedimentH;
		xsect->m_aFull = xsect->m_yFull * xsect->m_wMax;////矩形断面的面积
		xsect->m_rFull = xsect->m_aFull / (2.0 * (xsect->m_yFull + xsect->m_wMax));//矩形断面满流时的水力半径
		xsect->m_sFull = xsect->m_aFull * pow(xsect->m_rFull, 2. / 3.);//矩形断面满流时A*R^(2/3)的值，A为满流时过流断面面积，R为满流时水力半径
		aMax = RECT_ALFMAX * xsect->m_aFull;//充满度0.97时的面积
		xsect->m_sMax = aMax * pow(xsect->rect_closed_getRofA(aMax), 2. / 3.);//充满度0.97时A*R^(2/3)的值，A、R分别为充满度0.97时过流断面面积与水力半径
		xsect->m_ywMax = xsect->m_yFull;//断面最宽处的宽度
		xsect->m_aBot = xsect->m_yBot*xsect->m_wMax;//淤积的面积
		//根据淤积深度更新管段上下游端点处的offset
		Link[linkIndex].m_offset1 += Link[linkIndex].m_sedimentDeltaH;
		Link[linkIndex].m_offset2 += Link[linkIndex].m_sedimentDeltaH;
		break;
	case RECT_OPEN:
	case DUMMY:
	case CIRCULAR:
	case FORCE_MAIN:
	case EGGSHAPED:
	case HORSESHOE:
	case GOTHIC:
	case CATENARY:
	case SEMIELLIPTICAL:
	case BASKETHANDLE:
	case SEMICIRCULAR:	
	case RECT_TRIANG:
	case RECT_ROUND:
	case MOD_BASKET:
	case TRAPEZOIDAL:
	case TRIANGULAR:
	case PARABOLIC:
	case POWERFUNC:
	case HORIZ_ELLIPSE:
	case VERT_ELLIPSE:
	case ARCH:
		break;
	}
}
#pragma endregion

