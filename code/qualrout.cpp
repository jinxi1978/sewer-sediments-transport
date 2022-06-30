//-----------------------------------------------------------------------------
//   qualrout.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//             04/02/15   (Build 5.1.008)
//             04/30/15   (Build 5.1.009)
//             08/05/15   (Build 5.1.010)
//   Author:   L. Rossman
//
//   Water quality routing functions.
//
//   Build 5.1.008:
//   - Pollutant mass lost to seepage flow added to mass balance totals.
//   - Pollutant concen. increased when evaporation occurs.
//
//   Build 5.1.009:
//   - Criterion for dry link/storage node changed to avoid concen. blowup.
//
//   Build 5.1.010:
//   - Entire module re-written to be more compact and easier to follow.
//   - Neglible depth limit replaced with a negligible volume limit.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"
using namespace SWMMCPP;
void    projectClass::qualrout_init()
//
//  Input:   none
//  Output:  none
//  Purpose: initializes water quality concentrations in all nodes and links.
//
{
    int     i, p, isWet;
    double  c;

    for (i = 0; i < Nobjects[NODE]; i++)
    {
        isWet = ( Node[i].m_newDepth > FUDGE );
        for (p = 0; p < Nobjects[POLLUT]; p++)
        {
            c = 0.0;
            if ( isWet ) c = Pollut[p].m_initConcen;
            Node[i].m_oldQual[p] = c;
            Node[i].m_newQual[p] = c;
        }
    }

    for (i = 0; i < Nobjects[LINK]; i++)
    {
        isWet = ( Link[i].m_newDepth > FUDGE );
        for (p = 0; p < Nobjects[POLLUT]; p++)
        {
            c = 0.0;
            if ( isWet ) c = Pollut[p].m_initConcen;
            Link[i].m_oldQual[p] = c;
            Link[i].m_newQual[p] = c;
        }
    }
}

//=============================================================================

void projectClass::qualrout_execute(double tStep)
//
//  Input:   tStep = routing time step (sec)
//  Output:  none
//  Purpose: routes water quality constituents through the drainage
//           network over the current time step.
//
{
    int    i, j;
    double qIn, vAvg;

    // --- find mass flow each link contributes to its downstream node
    for ( i = 0; i < Nobjects[LINK]; i++ ) Link[i].findLinkMassFlow(Nobjects[POLLUT], tStep);

    // --- find new water quality concentration at each node  
    for (j = 0; j < Nobjects[NODE]; j++)
    {
        // --- get node inflow and average volume
        qIn = Node[j].m_inflow;
        vAvg = (Node[j].m_oldVolume + Node[j].m_newVolume) / 2.0;
        
        // --- save inflow concentrations if treatment applied
        if ( Node[j].m_treatment )
        {
            if ( qIn < ZERO ) qIn = 0.0;
            treatmnt_setInflow(qIn, Node[j].m_newQual);
        }
       
        // --- find new quality at the node 
        if ( Node[j].m_type == STORAGE || Node[j].m_oldVolume > FUDGE )
        {
            findStorageQual(j, tStep);
        }
        else Node[j].findNodeQual(Nobjects[POLLUT]);

        // --- apply treatment to new quality values
        if ( Node[j].m_treatment ) 
			treatmnt_treat(j, qIn, vAvg, tStep);
    }

    // --- find new water quality in each link
    for ( i = 0; i < Nobjects[LINK]; i++ ) findLinkQual(i, tStep);
}

//=============================================================================

double getMixedQual(double c, double v1, double wIn, double qIn, double tStep)
//
//  Input:   c = concentration in reactor at start of time step (mass/ft3)
//           v1 = volume in reactor at start of time step (ft3)
//           wIn = mass inflow rate (mass/sec)
//           qIn = flow inflow rate (cfs)
//           tStep = time step (sec)
//  Output:  returns pollutant concentration at end of time step (mass/ft3)
//  Purpose: finds pollutant concentration within a completely mixed reactor.
//
{
    double vIn, cIn, cMax;

    // --- if no inflow then reactor concentration is unchanged
    if ( qIn <= ZERO ) return c;

    // --- compute concentration of any inflow
    vIn = qIn * tStep;
    cIn = wIn * tStep / vIn;

    // --- mixture concen. can't exceed either original or inflow concen.
    cMax = MAX(c, cIn);

    // --- mix inflow with current reactor contents
    c = (c*v1 + wIn*tStep) / (v1 + vIn);
    c = MIN(c, cMax);
    c = MAX(c, 0.0);
    return c;
}


//=============================================================================

void linkClass::findLinkMassFlow(int pollutCount, double tStep)
//
//  Input:   i = link index
//           tStep = time step (sec)
//  Output:  none
//  Purpose: adds constituent mass flow out of link to the total
//           accumulation at the link's downstream node.
//
//  Note:    Node[].newQual[], the accumulator variable, already contains
//           contributions from runoff and other external inflows from
//           calculations made in routing_execute().
{
    int    p;
    double qLink, w;
	nodeClass* nodePtr;
    // --- find inflow to downstream node
    qLink = m_newFlow;

    // --- identify index of downstream node
	nodePtr = m_node2Ptr;
    if ( qLink < 0.0 ) nodePtr = m_node1Ptr;
    qLink = fabs(qLink);

    // --- examine each pollutant
    for (p = 0; p < pollutCount; p++)
    {
        // --- temporarily accumulate inflow load in Node[j].newQual
        w = qLink * m_oldQual[p];
        nodePtr->m_newQual[p] += w;

        // --- update total load transported by link
        m_totalLoad[p] += w * tStep;
    }
}

//=============================================================================

void nodeClass::findNodeQual(int pollutCount)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: finds new quality in a node with no storage volume.
//
{
    int    p;
    double qNode;

    // --- if there is flow into node then concen. = mass inflow/node flow
    qNode = m_inflow;
    if ( qNode > ZERO )
    {
        for (p = 0; p < pollutCount; p++)
        {
            m_newQual[p] /= qNode;
        }
    }

    // --- otherwise concen. is 0
    else for (p = 0; p < pollutCount; p++) m_newQual[p] = 0.0;
}

//=============================================================================

void projectClass::findLinkQual(int i, double tStep)
//
//  Input:   i = link index
//           tStep = routing time step (sec)
//  Output:  none
//  Purpose: finds new quality in a link at end of the current time step.
//
{
    int    j,                // upstream node index
           k,                // conduit index
           p;                // pollutant index
    double wIn,              // pollutant mass inflow rate (mass/sec)
           qIn,              // inflow rate (cfs)
           qSeep,            // rate of seepage loss (cfs)
           v1,               // link volume at start of time step (ft3)
           v2,               // link volume at end of time step (ft3)
           c1,               // current concentration within link (mass/ft3)
           c2,               // new concentration within link (mass/ft3)
           vEvap,            // volume lost to evaporation (ft3)
           vLosses,          // evap. + seepage volume loss (ft3)
           fEvap,            // evaporation concentration factor
           barrels;          // number of barrels in conduit

    // --- identify index of upstream node
    j = Link[i].m_node1;
    if ( Link[i].m_newFlow < 0.0 ) j = Link[i].m_node2;

    // --- link quality is that of upstream node when
    //     link is not a conduit or is a dummy link
    if ( Link[i].m_type != CONDUIT || Link[i].m_xsect.m_type == DUMMY )
    {
        for (p = 0; p < Nobjects[POLLUT]; p++)
        {
            Link[i].m_newQual[p] = Node[j].m_newQual[p];
        }
        return;
    }

    // --- get flow rates and evaporation loss
    k = Link[i].m_subIndex;
    barrels = Conduit[k].m_barrels;
    qIn  = fabs(Conduit[k].m_q1) * barrels;
    qSeep = Conduit[k].m_seepLossRate * barrels;
    vEvap = Conduit[k].m_evapLossRate * barrels * tStep;

    // --- get starting and ending volumes
    v1 = Link[i].m_oldVolume;
    v2 = Link[i].m_newVolume;
    vLosses = qSeep*tStep + vEvap;

    // --- compute factor by which concentrations are increased due to
    //     evaporation loss 
    fEvap = 1.0;
    if ( vEvap > 0.0 && v1 > ZeroVolume ) fEvap += vEvap / v1;

    // --- Steady Flow routing requires special treatment
    if ( RouteModel == SF )
    {
        findSFLinkQual(i, qSeep, fEvap, tStep);
        return;
    }

    // --- adjust inflow to compensate for volume change under Dynamic
    //     Wave routing (which produces just a single (out)flow rate
    //     for a conduit)
    if ( RouteModel == DW )
    {
        qIn = qIn + (v2 + vLosses - v1) / tStep; 
        qIn = MAX(qIn, 0.0);
    }

    // --- examine each pollutant
    for (p = 0; p < Nobjects[POLLUT]; p++)
    {
        // --- start with concen. at start of time step
        c1 = Link[i].m_oldQual[p];

        // --- update mass balance accounting for seepage loss
        massbal_addSeepageLoss(p, qSeep*c1);

        // --- increase concen. by evaporation factor
        c1 *= fEvap;

        // --- reduce concen. by 1st-order reaction
        c2 = getReactedQual(p, c1, v1, tStep);

        // --- mix resulting contents with inflow from upstream node
        wIn = Node[j].m_newQual[p]*qIn;
        c2 = getMixedQual(c2, v1, wIn, qIn, tStep);

        // --- set concen. to zero if remaining volume is negligible
        if ( v2 < ZeroVolume )
        {
            massbal_addToFinalStorage(p, c2 * v2);
            c2 = 0.0;
        }

        // --- assign new concen. to link
        Link[i].m_newQual[p] = c2;
    }
}

//=============================================================================

void  projectClass::findSFLinkQual(int i, double qSeep, double fEvap, double tStep)
//
//  Input:   i = link index
//           tStep = routing time step (sec)
//  Output:  none
//  Purpose: finds new quality in a link at end of the current time step for
//           Steady Flow routing.
//
{
    int j = Link[i].m_node1;
    int p;
    double c1, c2;
    double lossRate;

    // --- examine each pollutant
    for (p = 0; p < Nobjects[POLLUT]; p++)
    {
        // --- conduit's quality equals upstream node quality
        c1 = Node[j].m_newQual[p];

        // --- update mass balance accounting for seepage loss
        massbal_addSeepageLoss(p, qSeep*c1);

        // --- increase concen. by evaporation factor
        c1 *= fEvap;

        // --- apply first-order decay over travel time
        c2 = c1;
        if ( Pollut[p].m_kDecay > 0.0 )
        {
            c2 = c1 * exp(-Pollut[p].m_kDecay * tStep);
            c2 = MAX(0.0, c2);
            lossRate = (c1 - c2) * Link[i].m_newFlow;
            massbal_addReactedMass(p, lossRate);
        }
		Link[i].m_newQual[p] = c2;
    }
}

//=============================================================================

void  projectClass::findStorageQual(int j, double tStep)
//
//  Input:   j = node index
//           tStep = routing time step (sec)
//  Output:  none
//  Purpose: finds new quality in a node with storage volume.
//  
{
    int    p,                // pollutant index
           k;                // storage unit index
    double qIn,              // inflow rate (cfs)
           wIn,              // pollutant mass inflow rate (mass)
           v1,               // volume at start of time step (ft3)
           c1,               // initial pollutant concentration (mass/ft3)
           c2,               // final pollutant concentration (mass/ft3)
           qExfil = 0.0,     // exfiltration rate from storage unit (cfs)
           vEvap = 0.0,      // evaporation loss from storage unit (ft3)
           fEvap = 1.0;      // evaporation concentration factor

    // --- get inflow rate & initial volume
    qIn = Node[j].m_inflow;
    v1 = Node[j].m_oldVolume;

    // -- for storage nodes
    if ( Node[j].m_type == STORAGE )
    {    
        // --- update hydraulic residence time
        //     (HRT can be used in treatment functions)
		Node[j].updateHRT(Node[j].m_oldVolume, qIn, tStep);

        // --- get exfiltration rate and evaporation loss
        k = Node[j].m_subIndex;
        qExfil = Storage[k].m_exfilLoss / tStep;
        vEvap = Storage[k].m_evapLoss;

        // --- compute factor by which concentrations are increased due to
        //     evaporation loss (avoiding huge factors as storage unit
        //     dries out completely)
        if ( vEvap > 0.0 && v1 > ZeroVolume ) fEvap += vEvap / v1;
    }

    // --- for each pollutant
    for (p = 0; p < Nobjects[POLLUT]; p++)
    {
        // --- start with concen. at start of time step 
        c1 = Node[j].m_oldQual[p];

        // --- update mass balance accounting for exfiltration loss
        massbal_addSeepageLoss(p, qExfil*c1);

        // --- increase concen. by evaporation factor
        c1 *= fEvap;

        // --- apply first order reaction only if no separate treatment function
        if ( Node[j].m_treatment == NULL ||
             Node[j].m_treatment[p].m_equation == NULL )
        {
            c1 = getReactedQual(p, c1, v1, tStep);
        }

        // --- mix resulting contents with inflow from all sources
        //     (temporarily accumulated in Node[j].m_newQual)
        wIn = Node[j].m_newQual[p];
        c2 = getMixedQual(c1, v1, wIn, qIn, tStep);

		// --- set concen. to zero if remaining volume & inflow is negligible          //(5.1.015)
        if ( Node[j].m_newVolume <= ZeroVolume  && qIn <= FLOW_TOL)                //(5.1.015)
        {
            massbal_addToFinalStorage(p, c2 * Node[j].m_newVolume);
            c2 = 0.0;
        }

        // --- assign new concen. to node
        Node[j].m_newQual[p] = c2;
    }
}

//=============================================================================

void nodeClass::updateHRT(double v, double q, double tStep)
//
//  Input:   j = node index
//           v = storage volume (ft3)
//           q = inflow rate (cfs)
//           tStep = time step (sec)
//  Output:  none
//  Purpose: updates hydraulic residence time (i.e., water age) at a 
//           storage node.
//
{
    double hrt = m_storagePtr->m_hrt;
    if ( v < ZERO ) hrt = 0.0;
    else hrt = (hrt + tStep) * v / (v + q*tStep);
	m_storagePtr->m_hrt = MAX(hrt, 0.0);
}

//=============================================================================

double projectClass::getReactedQual(int p, double c, double v1, double tStep)
//
//  Input:   p = pollutant index
//           c = initial concentration (mass/ft3)
//           v1 = initial volume (ft3)
//           tStep = time step (sec)
//  Output:  none
//  Purpose: applies a first order reaction to a pollutant over a given
//           time step.
//
{
    double c2, lossRate;
    double kDecay = Pollut[p].m_kDecay;

    if ( kDecay == 0.0 ) return c;
    c2 = c * (1.0 - kDecay * tStep);
    c2 = MAX(0.0, c2);
    lossRate = (c - c2) * v1 / tStep;
    massbal_addReactedMass(p, lossRate);
    return c2;
}
#pragma region of jinxi
void linkClass::qualrout_sediment_circularFilled(int FlowUnits, int N_k_m, int sedimentPollutIndex, double sedimentDensity, double sedimentViodRatio,
	double alpha_sediment, double alpha_scour, double tStep, double sedimentV)
//根据管道流速及沉淀颗粒性质计算circularFilled断面的管段的淤积情况
//tStep表示计算的时间步长，s
//sedimentV表示颗粒的沉速，m/s
//i表示管段的索引号
{
	double depth, length, flow, velocity, volume, flowArea, core, hydR, k, m, Sstar, S, B, detaA, sedimentM, orgSedimentH, sedimentMax, sedimentAvaliable;
	xsectClass tempSect;
	//部分数值初始化
	m_sedimentDeltaH = 0;
	//将filled_circular断面的原始圆断面复制出来用来计算淤积后的淤积层高度,因为swmm现有的filled_circular对应
	//的断面计算函数精度有限，在淤积深度变化较小的情况下无法算出缩小淤积断面后对应的淤积高度，所以需要借用circular的断面计算函数
	tempSect.m_yFull = m_xsect.m_wMax;
	tempSect.m_aFull = PI / 4.0 * tempSect.m_yFull * tempSect.m_yFull;
	tempSect.m_type = CIRCULAR;
	//计算管段的平均输沙能力kg/m3
	//获得原始数据
	depth = m_newDepth*UCF(LENGTH);
	length = m_conduitPtr->m_length*UCF(LENGTH);//单位为m
	volume = m_newVolume*UCF(LENGTH)*UCF(LENGTH)*UCF(LENGTH);//将ft3转化为m3
	velocity = link_getVelocity(fabs(m_newFlow), m_newDepth)*UCF(LENGTH);//单位为m/s
	flow = fabs(m_newFlow * UCF(FLOW));//将流量单位转化为公制单位
	if (FlowUnits == 4) flow = flow*0.001;//LPS->MPS
	if (FlowUnits == 5) flow = flow*0.011574;//如果原始流量单位是MLD则转化MPS，后面计算中的流量单位均为MPS
	flowArea = m_xsect.xsect_getAofY(m_newDepth)*UCF(LENGTH)*UCF(LENGTH);//ft2转换为m2
	hydR = m_xsect.xsect_getRofY(m_newDepth)*UCF(LENGTH);//水力半径，ft转m
	orgSedimentH = m_sedimentH;//单位为ft
									 //计算V^3/gRω
	core = pow(velocity, 3) / (SI_GRAVITY*hydR*sedimentV);
	//calculate k and m by relationships between V^3/gRω and k、m
	k = m_xsect.linkScour_get_k(core, N_k_m);
	m = m_xsect.linkScour_get_m(core, N_k_m);
	//计算平均输沙能力kg/m3
	Sstar = k*pow(core, m);
	//计算管道平均含沙量kg/m3
	S = m_newQual[sedimentPollutIndex] * 0.001;//进行单位转换mg/L转换为kg/m3：newQual*0.001
	if (S >= Sstar)//沉积
	{
		//计算沉积断面平均宽度，即过水断面的平均宽度
		B = flowArea / depth;
		//计算淤积断面面积变化量
		detaA = B*alpha_sediment*sedimentV*(S - Sstar)*tStep / (sedimentDensity*(1 - sedimentViodRatio));//detaA单位为m2
		sedimentM = (1 - sedimentViodRatio)*sedimentDensity*detaA*length;//按照detaA需要沉淀的泥沙的质量kg
		sedimentMax = (S - Sstar)*volume;//允许沉积的最大质量
										 //计算沉淀以后管道内的污染物浓度，沉淀后的浓度不能小于管段的平均输沙能力，否则计算结果会产生震荡
										 //如果沉淀后污染物浓度小于管段平均携沙能力则将管段平均携沙能力作为沉积后的污染物浓度
		if (sedimentM > sedimentMax)//Sstar的单位由kg/m3转换为mg/L
		{
			//按照S-Sstar（即sedimentMax）计算需要沉淀的泥沙的质量kg
			//根据沉淀泥沙的质量求出沉淀泥沙导致的沉积物横断面面积变化
			detaA = sedimentMax / ((1 - sedimentViodRatio)*sedimentDensity*length);
			//将平均携沙能力作为管段的污染物浓度
			m_newQual[sedimentPollutIndex] = Sstar*1000.0;
			//根据上面的沉积物横断面面积的变化更新管段沉积物横断面面积、沉积物横断面高度及高度变化的结果
			m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
			m_sedimentH = tempSect.xsect_getYofA(m_sedimentArea);//单位为ft
			m_sedimentDeltaH = m_sedimentH - orgSedimentH;
		}
		//否则说明沉积后管段污染物浓度大于管段平均输沙能力，则按照沉积后管段污染物浓度计算detaA（上面已经算过了）及淤积深度的计算
		else
		{
			//只有新增的沉积物断面面积小于管道当前过流断面面积才进行沉积计算
			if (detaA / pow(UCF(LENGTH), 2) < m_xsect.m_aFull)
			{
				m_newQual[sedimentPollutIndex] = m_newQual[sedimentPollutIndex] - 1000 * sedimentM / volume;//单位mg/L;
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);//ft2
				m_sedimentH = tempSect.xsect_getYofA(m_sedimentArea);//单位为ft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;//单位为ft
			}
		}
	}
	else//冲刷S < Sstar
	{
		//冲刷时B应该为淤积断面顶宽m_xsect.sBot为filled_circular断面淤积部分顶层的宽度
		B = m_xsect.m_sBot*UCF(LENGTH);//单位转换为m
		detaA = B*alpha_scour*sedimentV*(S - Sstar)*tStep / (sedimentDensity*(1 - sedimentViodRatio));//detaA为负值
		sedimentM = (sedimentViodRatio - 1)*sedimentDensity*detaA*length;//按照detaA计算获得的冲刷带走的沉积物质量kg，是正值
		sedimentMax = (Sstar - S)*volume;//允许冲刷走的沉积物最大质量，冲刷后的管段污染物浓度不能超过水流的携沙能力Sstar，是正值
		sedimentAvaliable = m_sedimentArea*pow(UCF(LENGTH), 2)*length*(1 - sedimentViodRatio)*sedimentDensity;//管底沉积物的质量kg
																													//冲刷走的沉积物断面面积小于现有沉积物断面面积时
		if (m_sedimentArea + detaA / pow(UCF(LENGTH), 2) >= 0)
		{
			if (sedimentM >= sedimentMax)
			{
				//按照允许冲刷走的沉积物最大质量进行计算
				detaA = sedimentMax / ((sedimentViodRatio - 1)*sedimentDensity*length);//是负值
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
				m_sedimentH = tempSect.xsect_getYofA(m_sedimentArea);//单位为ft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;
				m_newQual[sedimentPollutIndex] = Sstar*1000.0;
			}
			else
			{
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
				m_sedimentH = tempSect.xsect_getYofA(m_sedimentArea);//单位为ft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;
				//计算泥沙冲刷后管道内污染物的浓度
				m_newQual[sedimentPollutIndex] = m_newQual[sedimentPollutIndex] + 1000 * sedimentM / volume;
			}
		}
		else if (m_sedimentArea > 0)//管底有沉积面积但是沉积面积小于冲刷造成的面积损失，则冲干净
		{
			if (sedimentAvaliable < sedimentMax)
			{
				//泥沙冲刷的质量kg
				sedimentM = (1 - sedimentViodRatio)*sedimentDensity*m_sedimentArea*length*pow(UCF(LENGTH), 2);
				//计算泥沙冲刷后管道内污染物的浓度
				m_newQual[sedimentPollutIndex] = m_newQual[sedimentPollutIndex] + 1000 * sedimentM / volume;
				m_sedimentArea = 0;
				m_sedimentH = 0;
				m_sedimentDeltaH = 0 - orgSedimentH;
			}
			else
			{
				//按照允许冲刷走的沉积物最大质量进行计算
				detaA = sedimentMax / ((sedimentViodRatio - 1)*sedimentDensity*length);
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
				m_sedimentH = tempSect.xsect_getYofA(m_sedimentArea);//单位为ft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;
				m_newQual[sedimentPollutIndex] = Sstar*1000.0;
			}
		}
	}
}

void linkClass::qualrout_sediment_rectClosedFilled(int FlowUnits,int N_k_m,int sedimentPollutIndex,double sedimentDensity,double sedimentViodRatio,
	double alpha_sediment,double alpha_scour,double tStep, double sedimentV)
//根据管道流速及沉淀颗粒性质计算circularFilled断面的管段的淤积情况
//tStep表示计算的时间步长，s
//sedimentV表示颗粒的沉速，m/s
//i表示管段的索引号
{
	double depth, length, flow, velocity, volume, flowArea, core, hydR, k, m, Sstar, S, B, detaA, sedimentM, orgSedimentH, sedimentMax, sedimentAvaliable;
	xsectClass tempSect;
	//部分数值初始化
	m_sedimentDeltaH = 0;
	//将断面的原始断面复制出来用来计算淤积后的淤积层高度
	tempSect.m_yFull = m_xsect.m_yFull + m_xsect.m_yBot;
	tempSect.m_aFull = tempSect.m_yFull * m_xsect.m_wMax;
	tempSect.m_type = RECT_CLOSED;
	//计算管段的平均输沙能力kg/m3
	//获得原始数据
	depth = m_newDepth*UCF(LENGTH);
	length = m_conduitPtr->m_length*UCF(LENGTH);//单位为m
	volume = m_newVolume*UCF(LENGTH)*UCF(LENGTH)*UCF(LENGTH);//将ft3转化为m3
	velocity = link_getVelocity(fabs(m_newFlow), m_newDepth)*UCF(LENGTH);//单位为m/s
	flow = fabs(m_newFlow * UCF(FLOW));//将流量单位转化为公制单位
	if (FlowUnits == 4) flow = flow*0.001;//LPS->MPS
	if (FlowUnits == 5) flow = flow*0.011574;//如果原始流量单位是MLD则转化MPS，后面计算中的流量单位均为MPS
	flowArea = m_xsect.xsect_getAofY(m_newDepth)*UCF(LENGTH)*UCF(LENGTH);//ft2转换为m2
	hydR = m_xsect.xsect_getRofY(m_newDepth)*UCF(LENGTH);//水力半径，ft转m
	orgSedimentH = m_sedimentH;//单位为ft
									 //计算V^3/gRω
	core = pow(velocity, 3) / (SI_GRAVITY*hydR*sedimentV);
	//根据V^3/gRω计算k、m
	k = m_xsect.linkScour_get_k(core, N_k_m);
	m = m_xsect.linkScour_get_m(core, N_k_m);
	//计算平均输沙能力kg/m3
	Sstar = k*pow(core, m);
	//计算管道平均含沙量kg/m3
	S = m_newQual[sedimentPollutIndex] * 0.001;//进行单位转换mg/L转换为kg/m3：newQual*0.001
	if (S >= Sstar)//沉积
	{
		//计算沉积断面平均宽度，即过水断面的平均宽度
		B = flowArea / depth;
		//计算淤积断面面积变化量
		detaA = B*alpha_sediment*sedimentV*(S - Sstar)*tStep / (sedimentDensity*(1 - sedimentViodRatio));//detaA单位为m2
		sedimentM = (1 - sedimentViodRatio)*sedimentDensity*detaA*length;//按照detaA需要沉淀的泥沙的质量kg
		sedimentMax = (S - Sstar)*volume;//允许沉积的最大质量
										 //计算沉淀以后管道内的污染物浓度，沉淀后的浓度不能小于管段的平均输沙能力，否则计算结果会产生震荡
										 //如果沉淀后污染物浓度小于管段平均携沙能力则将管段平均携沙能力作为沉积后的污染物浓度
		if (sedimentM > sedimentMax)//Sstar的单位由kg/m3转换为mg/L
		{
			//按照S-Sstar（即sedimentMax）计算需要沉淀的泥沙的质量kg
			//根据沉淀泥沙的质量求出沉淀泥沙导致的沉积物横断面面积变化
			detaA = sedimentMax / ((1 - sedimentViodRatio)*sedimentDensity*length);
			//将平均携沙能力作为管段的污染物浓度
			m_newQual[sedimentPollutIndex] = Sstar*1000.0;
			//根据上面的沉积物横断面面积的变化更新管段沉积物横断面面积、沉积物横断面高度及高度变化的结果
			m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
			m_sedimentH = m_sedimentArea/ m_xsect.m_wMax;//单位为ft
			m_sedimentDeltaH = m_sedimentH - orgSedimentH;
		}
		//否则说明沉积后管段污染物浓度大于管段平均输沙能力，则按照沉积后管段污染物浓度计算detaA（上面已经算过了）及淤积深度的计算
		else
		{
			//只有新增的沉积物断面面积小于管道当前过流断面面积才进行沉积计算
			if (detaA / pow(UCF(LENGTH), 2) < m_xsect.m_aFull)
			{
				m_newQual[sedimentPollutIndex] = m_newQual[sedimentPollutIndex] - 1000 * sedimentM / volume;//单位mg/L;
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);//ft2
				m_sedimentH = m_sedimentArea / m_xsect.m_wMax;//单位为ft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;//单位为ft
			}
		}
	}
	else//冲刷S < Sstar
	{
		//冲刷时B应该为淤积断面顶宽m_xsect.sBot为filled_circular断面淤积部分顶层的宽度
		B = m_xsect.m_wMax*UCF(LENGTH);//单位转换为m
		detaA = B*alpha_scour*sedimentV*(S - Sstar)*tStep / (sedimentDensity*(1 - sedimentViodRatio));//detaA为负值
		sedimentM = (sedimentViodRatio - 1)*sedimentDensity*detaA*length;//按照detaA计算获得的冲刷带走的沉积物质量kg，是正值
		sedimentMax = (Sstar - S)*volume;//允许冲刷走的沉积物最大质量，冲刷后的管段污染物浓度不能超过水流的携沙能力Sstar，是正值
		sedimentAvaliable = m_sedimentArea*pow(UCF(LENGTH), 2)*length*(1 - sedimentViodRatio)*sedimentDensity;//管底沉积物的质量kg
																													//冲刷走的沉积物断面面积小于现有沉积物断面面积时
		if (m_sedimentArea + detaA / pow(UCF(LENGTH), 2) >= 0)
		{
			if (sedimentM >= sedimentMax)
			{
				//按照允许冲刷走的沉积物最大质量进行计算
				detaA = sedimentMax / ((sedimentViodRatio - 1)*sedimentDensity*length);//是负值
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
				m_sedimentH = m_sedimentArea / m_xsect.m_wMax;//单位为ft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;
				m_newQual[sedimentPollutIndex] = Sstar*1000.0;
			}
			else
			{
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
				m_sedimentH = m_sedimentArea / m_xsect.m_wMax;//单位为ft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;
				//计算泥沙冲刷后管道内污染物的浓度
				m_newQual[sedimentPollutIndex] = m_newQual[sedimentPollutIndex] + 1000 * sedimentM / volume;
			}
		}
		else if (m_sedimentArea > 0)//管底有沉积面积但是沉积面积小于冲刷造成的面积损失，但是这时冲刷引起的面积变化没有考虑最大允许冲刷量，所以
			//还需要进一步判断
		{
			if (sedimentAvaliable < sedimentMax)//管底残留污染物量小于最大允许冲刷量（根据水流携沙能力和当前断面平均浓度求得最大允许冲刷量），则冲干净
			{
				//泥沙冲刷的质量kg
				sedimentM = (1 - sedimentViodRatio)*sedimentDensity*m_sedimentArea*length*pow(UCF(LENGTH), 2);
				//计算泥沙冲刷后管道内污染物的浓度
				m_newQual[sedimentPollutIndex] = m_newQual[sedimentPollutIndex] + 1000 * sedimentM / volume;
				m_sedimentArea = 0;
				m_sedimentH = 0;
				m_sedimentDeltaH = 0 - orgSedimentH;
			}
			else//管底残留污染物量大于最大允许冲刷量（根据水流携沙能力和当前断面平均浓度求得最大允许冲刷量），则会有残余管底沉积物
			{
				//按照允许冲刷走的沉积物最大质量进行计算
				detaA = sedimentMax / ((sedimentViodRatio - 1)*sedimentDensity*length);
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
				m_sedimentH = m_sedimentArea / m_xsect.m_wMax;//单位为ft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;
				m_newQual[sedimentPollutIndex] = Sstar*1000.0;
			}
		}
	}
}

//this function is called in routing_execute( ) like this	
//void projectClass::routing_execute(int routingModel, double routingStep)
//{................
//        // --- route quality through the drainage network
//        if ( Nobjects[POLLUT] > 0 && !IgnoreQuality ) 
//        {
//            qualrout_execute(routingStep);
//        }
//		//added by jinxi进行管段沉积物淤积冲刷的模拟
//		if (Nobjects[POLLUT] > 0 && sedimentPollutIndex < Nobjects[POLLUT] && !IgnoreSediment)
//		{
//			qualrout_sediment(routingStep);
//			for (j = 0; j < Nobjects[LINK]; j++)
//				xsect_updateCrossSection(j);
//		}
//		//end of content added by jinxi 
//        // --- remove evaporation, infiltration & outflows from system
//        removeStorageLosses(routingStep);
//        removeConduitLosses();
//       removeOutflows(routingStep);
//................}
void	projectClass::qualrout_sediment(double tStep)
{
	int i=0;
	double sedimentV;
	//首先计算颗粒沉速，通过之前输入的粒径与沉速之间的关系通过插值法求得沉速
	if (sedimentPollutIndex<0 || sedimentPollutIndex > Nobjects[POLLUT]) return;
	sedimentV = Link[i].m_xsect.linkScour_get_V(sedimentParticalDia, N_D_V);
	for (i = 0; i < Nobjects[LINK]; i++)
	{
		if (Link[i].m_type == CONDUIT && Link[i].m_xsect.m_type == FILLED_CIRCULAR)
		{
			Link[i].qualrout_sediment_circularFilled(tStep, N_k_m, sedimentPollutIndex, sedimentDensity, sedimentViodRatio,
				alpha_sediment, alpha_scour, tStep, sedimentV);
		}
		else if (Link[i].m_type == CONDUIT && Link[i].m_xsect.m_type == RECT_CLOSED)
		{
			Link[i].qualrout_sediment_rectClosedFilled(tStep,N_k_m,sedimentPollutIndex,sedimentDensity,sedimentViodRatio,
				alpha_sediment,alpha_scour,tStep,sedimentV);
		}
	}
}


#pragma endregion

