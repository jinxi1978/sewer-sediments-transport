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
//���ݹܵ����ټ�����������ʼ���circularFilled����Ĺܶε��ٻ����
//tStep��ʾ�����ʱ�䲽����s
//sedimentV��ʾ�����ĳ��٣�m/s
//i��ʾ�ܶε�������
{
	double depth, length, flow, velocity, volume, flowArea, core, hydR, k, m, Sstar, S, B, detaA, sedimentM, orgSedimentH, sedimentMax, sedimentAvaliable;
	xsectClass tempSect;
	//������ֵ��ʼ��
	m_sedimentDeltaH = 0;
	//��filled_circular�����ԭʼԲ���渴�Ƴ������������ٻ�����ٻ���߶�,��Ϊswmm���е�filled_circular��Ӧ
	//�Ķ�����㺯���������ޣ����ٻ���ȱ仯��С��������޷������С�ٻ�������Ӧ���ٻ��߶ȣ�������Ҫ����circular�Ķ�����㺯��
	tempSect.m_yFull = m_xsect.m_wMax;
	tempSect.m_aFull = PI / 4.0 * tempSect.m_yFull * tempSect.m_yFull;
	tempSect.m_type = CIRCULAR;
	//����ܶε�ƽ����ɳ����kg/m3
	//���ԭʼ����
	depth = m_newDepth*UCF(LENGTH);
	length = m_conduitPtr->m_length*UCF(LENGTH);//��λΪm
	volume = m_newVolume*UCF(LENGTH)*UCF(LENGTH)*UCF(LENGTH);//��ft3ת��Ϊm3
	velocity = link_getVelocity(fabs(m_newFlow), m_newDepth)*UCF(LENGTH);//��λΪm/s
	flow = fabs(m_newFlow * UCF(FLOW));//��������λת��Ϊ���Ƶ�λ
	if (FlowUnits == 4) flow = flow*0.001;//LPS->MPS
	if (FlowUnits == 5) flow = flow*0.011574;//���ԭʼ������λ��MLD��ת��MPS����������е�������λ��ΪMPS
	flowArea = m_xsect.xsect_getAofY(m_newDepth)*UCF(LENGTH)*UCF(LENGTH);//ft2ת��Ϊm2
	hydR = m_xsect.xsect_getRofY(m_newDepth)*UCF(LENGTH);//ˮ���뾶��ftתm
	orgSedimentH = m_sedimentH;//��λΪft
									 //����V^3/gR��
	core = pow(velocity, 3) / (SI_GRAVITY*hydR*sedimentV);
	//calculate k and m by relationships between V^3/gR�� and k��m
	k = m_xsect.linkScour_get_k(core, N_k_m);
	m = m_xsect.linkScour_get_m(core, N_k_m);
	//����ƽ����ɳ����kg/m3
	Sstar = k*pow(core, m);
	//����ܵ�ƽ����ɳ��kg/m3
	S = m_newQual[sedimentPollutIndex] * 0.001;//���е�λת��mg/Lת��Ϊkg/m3��newQual*0.001
	if (S >= Sstar)//����
	{
		//�����������ƽ����ȣ�����ˮ�����ƽ�����
		B = flowArea / depth;
		//�����ٻ���������仯��
		detaA = B*alpha_sediment*sedimentV*(S - Sstar)*tStep / (sedimentDensity*(1 - sedimentViodRatio));//detaA��λΪm2
		sedimentM = (1 - sedimentViodRatio)*sedimentDensity*detaA*length;//����detaA��Ҫ�������ɳ������kg
		sedimentMax = (S - Sstar)*volume;//����������������
										 //��������Ժ�ܵ��ڵ���Ⱦ��Ũ�ȣ�������Ũ�Ȳ���С�ڹܶε�ƽ����ɳ����������������������
										 //����������Ⱦ��Ũ��С�ڹܶ�ƽ��Яɳ�����򽫹ܶ�ƽ��Яɳ������Ϊ���������Ⱦ��Ũ��
		if (sedimentM > sedimentMax)//Sstar�ĵ�λ��kg/m3ת��Ϊmg/L
		{
			//����S-Sstar����sedimentMax��������Ҫ�������ɳ������kg
			//���ݳ�����ɳ���������������ɳ���µĳ�������������仯
			detaA = sedimentMax / ((1 - sedimentViodRatio)*sedimentDensity*length);
			//��ƽ��Яɳ������Ϊ�ܶε���Ⱦ��Ũ��
			m_newQual[sedimentPollutIndex] = Sstar*1000.0;
			//��������ĳ�������������ı仯���¹ܶγ������������������������߶ȼ��߶ȱ仯�Ľ��
			m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
			m_sedimentH = tempSect.xsect_getYofA(m_sedimentArea);//��λΪft
			m_sedimentDeltaH = m_sedimentH - orgSedimentH;
		}
		//����˵��������ܶ���Ⱦ��Ũ�ȴ��ڹܶ�ƽ����ɳ���������ճ�����ܶ���Ⱦ��Ũ�ȼ���detaA�������Ѿ�����ˣ����ٻ���ȵļ���
		else
		{
			//ֻ�������ĳ�����������С�ڹܵ���ǰ������������Ž��г�������
			if (detaA / pow(UCF(LENGTH), 2) < m_xsect.m_aFull)
			{
				m_newQual[sedimentPollutIndex] = m_newQual[sedimentPollutIndex] - 1000 * sedimentM / volume;//��λmg/L;
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);//ft2
				m_sedimentH = tempSect.xsect_getYofA(m_sedimentArea);//��λΪft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;//��λΪft
			}
		}
	}
	else//��ˢS < Sstar
	{
		//��ˢʱBӦ��Ϊ�ٻ����涥��m_xsect.sBotΪfilled_circular�����ٻ����ֶ���Ŀ��
		B = m_xsect.m_sBot*UCF(LENGTH);//��λת��Ϊm
		detaA = B*alpha_scour*sedimentV*(S - Sstar)*tStep / (sedimentDensity*(1 - sedimentViodRatio));//detaAΪ��ֵ
		sedimentM = (sedimentViodRatio - 1)*sedimentDensity*detaA*length;//����detaA�����õĳ�ˢ���ߵĳ���������kg������ֵ
		sedimentMax = (Sstar - S)*volume;//�����ˢ�ߵĳ����������������ˢ��Ĺܶ���Ⱦ��Ũ�Ȳ��ܳ���ˮ����Яɳ����Sstar������ֵ
		sedimentAvaliable = m_sedimentArea*pow(UCF(LENGTH), 2)*length*(1 - sedimentViodRatio)*sedimentDensity;//�ܵ׳����������kg
																													//��ˢ�ߵĳ�����������С�����г�����������ʱ
		if (m_sedimentArea + detaA / pow(UCF(LENGTH), 2) >= 0)
		{
			if (sedimentM >= sedimentMax)
			{
				//���������ˢ�ߵĳ���������������м���
				detaA = sedimentMax / ((sedimentViodRatio - 1)*sedimentDensity*length);//�Ǹ�ֵ
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
				m_sedimentH = tempSect.xsect_getYofA(m_sedimentArea);//��λΪft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;
				m_newQual[sedimentPollutIndex] = Sstar*1000.0;
			}
			else
			{
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
				m_sedimentH = tempSect.xsect_getYofA(m_sedimentArea);//��λΪft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;
				//������ɳ��ˢ��ܵ�����Ⱦ���Ũ��
				m_newQual[sedimentPollutIndex] = m_newQual[sedimentPollutIndex] + 1000 * sedimentM / volume;
			}
		}
		else if (m_sedimentArea > 0)//�ܵ��г���������ǳ������С�ڳ�ˢ��ɵ������ʧ�����ɾ�
		{
			if (sedimentAvaliable < sedimentMax)
			{
				//��ɳ��ˢ������kg
				sedimentM = (1 - sedimentViodRatio)*sedimentDensity*m_sedimentArea*length*pow(UCF(LENGTH), 2);
				//������ɳ��ˢ��ܵ�����Ⱦ���Ũ��
				m_newQual[sedimentPollutIndex] = m_newQual[sedimentPollutIndex] + 1000 * sedimentM / volume;
				m_sedimentArea = 0;
				m_sedimentH = 0;
				m_sedimentDeltaH = 0 - orgSedimentH;
			}
			else
			{
				//���������ˢ�ߵĳ���������������м���
				detaA = sedimentMax / ((sedimentViodRatio - 1)*sedimentDensity*length);
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
				m_sedimentH = tempSect.xsect_getYofA(m_sedimentArea);//��λΪft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;
				m_newQual[sedimentPollutIndex] = Sstar*1000.0;
			}
		}
	}
}

void linkClass::qualrout_sediment_rectClosedFilled(int FlowUnits,int N_k_m,int sedimentPollutIndex,double sedimentDensity,double sedimentViodRatio,
	double alpha_sediment,double alpha_scour,double tStep, double sedimentV)
//���ݹܵ����ټ�����������ʼ���circularFilled����Ĺܶε��ٻ����
//tStep��ʾ�����ʱ�䲽����s
//sedimentV��ʾ�����ĳ��٣�m/s
//i��ʾ�ܶε�������
{
	double depth, length, flow, velocity, volume, flowArea, core, hydR, k, m, Sstar, S, B, detaA, sedimentM, orgSedimentH, sedimentMax, sedimentAvaliable;
	xsectClass tempSect;
	//������ֵ��ʼ��
	m_sedimentDeltaH = 0;
	//�������ԭʼ���渴�Ƴ������������ٻ�����ٻ���߶�
	tempSect.m_yFull = m_xsect.m_yFull + m_xsect.m_yBot;
	tempSect.m_aFull = tempSect.m_yFull * m_xsect.m_wMax;
	tempSect.m_type = RECT_CLOSED;
	//����ܶε�ƽ����ɳ����kg/m3
	//���ԭʼ����
	depth = m_newDepth*UCF(LENGTH);
	length = m_conduitPtr->m_length*UCF(LENGTH);//��λΪm
	volume = m_newVolume*UCF(LENGTH)*UCF(LENGTH)*UCF(LENGTH);//��ft3ת��Ϊm3
	velocity = link_getVelocity(fabs(m_newFlow), m_newDepth)*UCF(LENGTH);//��λΪm/s
	flow = fabs(m_newFlow * UCF(FLOW));//��������λת��Ϊ���Ƶ�λ
	if (FlowUnits == 4) flow = flow*0.001;//LPS->MPS
	if (FlowUnits == 5) flow = flow*0.011574;//���ԭʼ������λ��MLD��ת��MPS����������е�������λ��ΪMPS
	flowArea = m_xsect.xsect_getAofY(m_newDepth)*UCF(LENGTH)*UCF(LENGTH);//ft2ת��Ϊm2
	hydR = m_xsect.xsect_getRofY(m_newDepth)*UCF(LENGTH);//ˮ���뾶��ftתm
	orgSedimentH = m_sedimentH;//��λΪft
									 //����V^3/gR��
	core = pow(velocity, 3) / (SI_GRAVITY*hydR*sedimentV);
	//����V^3/gR�ؼ���k��m
	k = m_xsect.linkScour_get_k(core, N_k_m);
	m = m_xsect.linkScour_get_m(core, N_k_m);
	//����ƽ����ɳ����kg/m3
	Sstar = k*pow(core, m);
	//����ܵ�ƽ����ɳ��kg/m3
	S = m_newQual[sedimentPollutIndex] * 0.001;//���е�λת��mg/Lת��Ϊkg/m3��newQual*0.001
	if (S >= Sstar)//����
	{
		//�����������ƽ����ȣ�����ˮ�����ƽ�����
		B = flowArea / depth;
		//�����ٻ���������仯��
		detaA = B*alpha_sediment*sedimentV*(S - Sstar)*tStep / (sedimentDensity*(1 - sedimentViodRatio));//detaA��λΪm2
		sedimentM = (1 - sedimentViodRatio)*sedimentDensity*detaA*length;//����detaA��Ҫ�������ɳ������kg
		sedimentMax = (S - Sstar)*volume;//����������������
										 //��������Ժ�ܵ��ڵ���Ⱦ��Ũ�ȣ�������Ũ�Ȳ���С�ڹܶε�ƽ����ɳ����������������������
										 //����������Ⱦ��Ũ��С�ڹܶ�ƽ��Яɳ�����򽫹ܶ�ƽ��Яɳ������Ϊ���������Ⱦ��Ũ��
		if (sedimentM > sedimentMax)//Sstar�ĵ�λ��kg/m3ת��Ϊmg/L
		{
			//����S-Sstar����sedimentMax��������Ҫ�������ɳ������kg
			//���ݳ�����ɳ���������������ɳ���µĳ�������������仯
			detaA = sedimentMax / ((1 - sedimentViodRatio)*sedimentDensity*length);
			//��ƽ��Яɳ������Ϊ�ܶε���Ⱦ��Ũ��
			m_newQual[sedimentPollutIndex] = Sstar*1000.0;
			//��������ĳ�������������ı仯���¹ܶγ������������������������߶ȼ��߶ȱ仯�Ľ��
			m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
			m_sedimentH = m_sedimentArea/ m_xsect.m_wMax;//��λΪft
			m_sedimentDeltaH = m_sedimentH - orgSedimentH;
		}
		//����˵��������ܶ���Ⱦ��Ũ�ȴ��ڹܶ�ƽ����ɳ���������ճ�����ܶ���Ⱦ��Ũ�ȼ���detaA�������Ѿ�����ˣ����ٻ���ȵļ���
		else
		{
			//ֻ�������ĳ�����������С�ڹܵ���ǰ������������Ž��г�������
			if (detaA / pow(UCF(LENGTH), 2) < m_xsect.m_aFull)
			{
				m_newQual[sedimentPollutIndex] = m_newQual[sedimentPollutIndex] - 1000 * sedimentM / volume;//��λmg/L;
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);//ft2
				m_sedimentH = m_sedimentArea / m_xsect.m_wMax;//��λΪft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;//��λΪft
			}
		}
	}
	else//��ˢS < Sstar
	{
		//��ˢʱBӦ��Ϊ�ٻ����涥��m_xsect.sBotΪfilled_circular�����ٻ����ֶ���Ŀ��
		B = m_xsect.m_wMax*UCF(LENGTH);//��λת��Ϊm
		detaA = B*alpha_scour*sedimentV*(S - Sstar)*tStep / (sedimentDensity*(1 - sedimentViodRatio));//detaAΪ��ֵ
		sedimentM = (sedimentViodRatio - 1)*sedimentDensity*detaA*length;//����detaA�����õĳ�ˢ���ߵĳ���������kg������ֵ
		sedimentMax = (Sstar - S)*volume;//�����ˢ�ߵĳ����������������ˢ��Ĺܶ���Ⱦ��Ũ�Ȳ��ܳ���ˮ����Яɳ����Sstar������ֵ
		sedimentAvaliable = m_sedimentArea*pow(UCF(LENGTH), 2)*length*(1 - sedimentViodRatio)*sedimentDensity;//�ܵ׳����������kg
																													//��ˢ�ߵĳ�����������С�����г�����������ʱ
		if (m_sedimentArea + detaA / pow(UCF(LENGTH), 2) >= 0)
		{
			if (sedimentM >= sedimentMax)
			{
				//���������ˢ�ߵĳ���������������м���
				detaA = sedimentMax / ((sedimentViodRatio - 1)*sedimentDensity*length);//�Ǹ�ֵ
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
				m_sedimentH = m_sedimentArea / m_xsect.m_wMax;//��λΪft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;
				m_newQual[sedimentPollutIndex] = Sstar*1000.0;
			}
			else
			{
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
				m_sedimentH = m_sedimentArea / m_xsect.m_wMax;//��λΪft
				m_sedimentDeltaH = m_sedimentH - orgSedimentH;
				//������ɳ��ˢ��ܵ�����Ⱦ���Ũ��
				m_newQual[sedimentPollutIndex] = m_newQual[sedimentPollutIndex] + 1000 * sedimentM / volume;
			}
		}
		else if (m_sedimentArea > 0)//�ܵ��г���������ǳ������С�ڳ�ˢ��ɵ������ʧ��������ʱ��ˢ���������仯û�п�����������ˢ��������
			//����Ҫ��һ���ж�
		{
			if (sedimentAvaliable < sedimentMax)//�ܵײ�����Ⱦ����С����������ˢ��������ˮ��Яɳ�����͵�ǰ����ƽ��Ũ�������������ˢ���������ɾ�
			{
				//��ɳ��ˢ������kg
				sedimentM = (1 - sedimentViodRatio)*sedimentDensity*m_sedimentArea*length*pow(UCF(LENGTH), 2);
				//������ɳ��ˢ��ܵ�����Ⱦ���Ũ��
				m_newQual[sedimentPollutIndex] = m_newQual[sedimentPollutIndex] + 1000 * sedimentM / volume;
				m_sedimentArea = 0;
				m_sedimentH = 0;
				m_sedimentDeltaH = 0 - orgSedimentH;
			}
			else//�ܵײ�����Ⱦ����������������ˢ��������ˮ��Яɳ�����͵�ǰ����ƽ��Ũ�������������ˢ����������в���ܵ׳�����
			{
				//���������ˢ�ߵĳ���������������м���
				detaA = sedimentMax / ((sedimentViodRatio - 1)*sedimentDensity*length);
				m_sedimentArea += detaA / pow(UCF(LENGTH), 2);
				m_sedimentH = m_sedimentArea / m_xsect.m_wMax;//��λΪft
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
//		//added by jinxi���йܶγ������ٻ���ˢ��ģ��
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
	//���ȼ���������٣�ͨ��֮ǰ��������������֮��Ĺ�ϵͨ����ֵ����ó���
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

