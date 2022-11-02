package multi_energy_dvcs_bsa;

import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;

import org.freehep.math.minuit.FCNBase;
import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMigrad;
import org.freehep.math.minuit.MnStrategy;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.utils.groups.IndexedList;

import org.jlab.io.evio.EvioDataBank;
import org.jlab.io.hipo.HipoDataBank;

public class p_theta_e_elastic_correction {
	
	static HipoDataSource reader = new HipoDataSource();
	static HipoDataSource readerc = new HipoDataSource();
	static HipoDataSource readerc1 = new HipoDataSource();
	
	static H1F hWs = new H1F("hWs", "hWs", 200, 0.7, 1.2);
	static H1F hpcorWs = new H1F("hpcorWs", "hpcorWs", 200, 0.7, 1.2);
	static H1F hthpcorWs = new H1F("hpcorWs1", "hpcorWs1", 200, 0.7, 1.2);
	
	static IndexedList<H1F> histGroups_W_sec = new IndexedList<H1F>(1);
	
	public static void histos_W_sec() {
		for(int ihist1 = 0; ihist1 < 6; ihist1++) {
			H1F hW_sec = new H1F("hW_sec", "hW_sec", 200, 0.7, 1.2);
			histGroups_W_sec.add(hW_sec, ihist1);
		}
	}
	
	static IndexedList<H1F> histGroups_p_cor_W_sec = new IndexedList<H1F>(1);
	
	public static void histos_p_cor_W_sec() {
		for(int ihist1cor = 0; ihist1cor < 6; ihist1cor++) {
			H1F hp_cor_W_sec = new H1F("hp_cor_W_sec", "hp_cor_W_sec", 200, 0.7, 1.2);
			histGroups_p_cor_W_sec.add(hp_cor_W_sec, ihist1cor);
		}
	}
	
	static IndexedList<H1F> histGroups_theta_p_cor_W_sec = new IndexedList<H1F>(1);
	
	public static void histos_theta_p_cor_W_sec() {
		for(int ihist1cor1 = 0; ihist1cor1 < 6; ihist1cor1++) {
			H1F htheta_p_cor_W_sec = new H1F("htheta_p_cor_W_sec", "htheta_p_cor_W_sec", 200, 0.7, 1.2);
			histGroups_theta_p_cor_W_sec.add(htheta_p_cor_W_sec, ihist1cor1);
		}
	}
	
	static IndexedList<H2F> histGroups_W_vs_phi = new IndexedList<H2F>(1);
	
	public static void histos2() {
		for(int ihist2 = 0; ihist2 < 9; ihist2++) {
			H2F hW_vs_phi = new H2F("hW_vs_phi", "hW_vs_phi", 600, -150, 210, 200, 0.7, 1.2);
			histGroups_W_vs_phi.add(hW_vs_phi, ihist2);
		}
	}
	
	static IndexedList<H1F> histGroups_W_theta_phi_bin = new IndexedList<H1F>(3);
	
	public static void histos_W_vs_phi_bin() {
		for(int ihistsector = 0; ihistsector < 6; ihistsector++) {
			for(int ihisttheta = 0; ihisttheta < 9; ihisttheta++) {
				for(int ihistphi = 0; ihistphi < 4; ihistphi++) {
					H1F hW_vs_phi_bin = new H1F("hW_vs_phi_bin", "hW_vs_phi_bin", 200, 0.7, 1.2);
					histGroups_W_theta_phi_bin.add(hW_vs_phi_bin, ihistsector, ihisttheta, ihistphi);
				}
			}
		}
	}
	
	static IndexedList<H1F> histGroups_p_cor_W_theta_phi_bin = new IndexedList<H1F>(3);
	
	public static void histos_p_cor_W_vs_phi_bin() {
		for(int ihistsectorc = 0; ihistsectorc < 6; ihistsectorc++) {
			for(int ihistthetac = 0; ihistthetac < 9; ihistthetac++) {
				for(int ihistphic = 0; ihistphic < 4; ihistphic++) {
					H1F hp_cor_W_theta_phi_bin = new H1F("hp_cor_W_theta_phi_bin", "hp_cor_W_theta_phi_bin", 200, 0.7, 1.2);
					histGroups_p_cor_W_theta_phi_bin.add(hp_cor_W_theta_phi_bin, ihistsectorc, ihistthetac, ihistphic);
				}
			}
		}
	}
	
	static IndexedList<H1F> histGroups_theta_p_cor_W_theta_phi_bin = new IndexedList<H1F>(3);
	
	public static void histos_theta_p_cor_W_theta_phi_bin() {
		for(int ihistsectorc1 = 0; ihistsectorc1 < 6; ihistsectorc1++) {
			for(int ihistthetac1 = 0; ihistthetac1 < 9; ihistthetac1++) {
				for(int ihistphic1 = 0; ihistphic1 < 4; ihistphic1++) {
					H1F htheta_p_cor_W_theta_phi_bin = new H1F("htheta_p_cor_W_theta_phi_bin", "htheta_p_cor_W_theta_phi_bin", 200, 0.7, 1.2);
					histGroups_theta_p_cor_W_theta_phi_bin.add(htheta_p_cor_W_theta_phi_bin, ihistsectorc1, ihistthetac1, ihistphic1);
				}
			}
		}
	}
	
	static IndexedList<H1F> histGroups_dp = new IndexedList<H1F>(3);
	
	public static void histos_dp() {
		for(int ihistsectordp = 0; ihistsectordp < 6; ihistsectordp++) {
			for(int ihistthetadp = 0; ihistthetadp < 9; ihistthetadp++) {
				for(int ihistphidp = 0; ihistphidp < 4; ihistphidp++) {
					H1F hdp = new H1F("hdp", "hdp", 200, 0.85, 1.15);
					histGroups_dp.add(hdp, ihistsectordp, ihistthetadp, ihistphidp);
				}
			}
		}
	}
	
	static IndexedList<H1F> histGroups_dtheta = new IndexedList<H1F>(3);
	
	public static void histos_dtheta() {
		for(int ihistsectordth = 0; ihistsectordth < 6; ihistsectordth++) {
			for(int ihistthetadth = 0; ihistthetadth < 9; ihistthetadth++) {
				for(int ihistphidth = 0; ihistphidth < 4; ihistphidth++) {
					H1F hdtheta = new H1F("hdtheta", "hdtheta", 200, 0, 2);
					histGroups_dtheta.add(hdtheta, ihistsectordth, ihistthetadth, ihistphidth);
				}
			}
		}
	}

	static void processEvent(DataEvent event, int eventCounter, byte setsector) {
		
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Track") && event.hasBank("REC::Traj"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int ecount = 0;
			int protoncount = 0;
			int e_pindex = -1;
			int proton_pindex = -1;
			float M = (float) 0.938;
			float E = (float) 6.535;
			float W = -20;
			float dp = -1;
			float theta_deg_e = 0;
			float phi_deg_e = -200;
			float phi_rot_deg_e = -200;
			double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0);
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0);
			Vector3D v_e = new Vector3D (0, 0, 0);
			Vector3D v_proton = new Vector3D (0, 0, 0);
			float cdchi2_proton = 100;
			short cdNDF_proton = 100;
			Vector3D dc_hit_e = new Vector3D (0, 0, 0);
			Vector3D dc_hit_rot_e = new Vector3D (0, 0, 0);
			for(int i = 0; i < rec.rows(); i++)
			{
				int pid = rec.getInt("pid", i);
				float px = rec.getFloat("px", i);
				float py = rec.getFloat("py", i);
				float pz = rec.getFloat("pz", i);
				float vx = rec.getFloat("vx", i);
				float vy = rec.getFloat("vy", i);
				float vz = rec.getFloat("vz", i);
				float p = (float) Math.sqrt((px*px)+(py*py)+(pz*pz));
				if(pid == 11)
				{
					ecount++;
					if(ecount == 1)
					{
						p_e = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.000511*0.000511)));
						v_e = new Vector3D (vx, vy, vz);
						W = (float) Math.sqrt((M*M)+(2*M*(E-p))+(2*E*p*((pz/p)-1)));
						dp = (float) (E/(p*(1+(2*E*(Math.sin(p_e.theta()/2))*(Math.sin(p_e.theta()/2))/M))));
						e_pindex = i;
					}
				}
				if(pid == 2212)
				{
					protoncount++;
					if(protoncount == 1)
					{
						p_proton = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.938272*0.938272)));
						v_proton = new Vector3D (vx, vy, vz);
						proton_pindex = i;
					}
				}
			}
			HipoDataBank rectrac = (HipoDataBank) event.getBank("REC::Track");
			byte e_detector = 0;
			byte e_sector = 0;
			byte proton_detector = 0;
			for(int j = 0; j < rectrac.rows(); j++)
			{
				short pindex = rectrac.getShort("pindex", j);
				byte detector = rectrac.getByte("detector", j);
				byte sector = rectrac.getByte("sector", j);
				float chi2 = rectrac.getFloat("chi2", j);
				short NDF = rectrac.getShort("NDF", j);
				if(pindex == e_pindex)
				{
					e_detector = detector;
					e_sector = sector;
				}
				if(pindex == proton_pindex)
				{
					proton_detector = detector;
					cdchi2_proton = chi2;
					cdNDF_proton = NDF;
				}
			}
			HipoDataBank rectraj = (HipoDataBank) event.getBank("REC::Traj");
			for(int k = 0; k < rectraj.rows(); k++)
			{
				short pindex = rectraj.getShort("pindex", k);
				byte detector = rectraj.getByte("detector", k);
				byte layer = rectraj.getByte("layer", k);
				float x = rectraj.getFloat("x", k);
				float y = rectraj.getFloat("y", k);
				float z = rectraj.getFloat("z", k);
				if(pindex == e_pindex && detector == 6 && e_detector == 6 && layer ==  6)
				{
					dc_hit_e = new Vector3D (x, y, z);
					theta_deg_e = (float) (57.3*dc_hit_e.theta());
					phi_deg_e = (float) (57.3*dc_hit_e.phi());
					if(phi_deg_e <= -150) phi_deg_e = phi_deg_e+360;
					dc_hit_rot_e = new Vector3D (x*Math.cos(phi_rot[e_sector-1]/57.3)+y*Math.sin(phi_rot[e_sector-1]/57.3),
							y*Math.cos(phi_rot[e_sector-1]/57.3)-x*Math.sin(phi_rot[e_sector-1]/57.3), z);
					phi_rot_deg_e = (float) (57.3*dc_hit_rot_e.phi());
				}
			}
			
			if(ecount == 1 && protoncount == 1
					&& v_e.z() > -12 && v_e.z() < 7 && Math.abs(v_e.z()-v_proton.z()) < 2.5 + (2.5/p_proton.p())
					&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30)))
				//	&& 57.3*p_proton.theta() < 75)
			{
				hWs.fill(W);
				histGroups_W_sec.getItem(e_sector-1).fill(W);
				
				double[] theta_bnd  = new double[] {7, 8.5, 10, 12, 14, 16, 18, 21, 24, 33};
				for(int thbin = 0; thbin < 9; thbin++)
				{
					if(theta_deg_e > theta_bnd[thbin] && theta_deg_e <= theta_bnd[thbin+1])
					{
						histGroups_W_vs_phi.getItem(thbin).fill(phi_deg_e, W);
						
						double[] phi_rot_bnd = new double[] {-30, -9, 0, 9, 30};
						for(int phbin = 0; phbin < 4; phbin++)
						{
							if(phi_rot_deg_e > phi_rot_bnd[phbin] && phi_rot_deg_e <= phi_rot_bnd[phbin+1])
							{
								histGroups_W_theta_phi_bin.getItem(e_sector-1, thbin, phbin).fill(W);
								histGroups_dp.getItem(e_sector-1, thbin, phbin).fill(dp);
								break;
							}
						}
						break;
					}
				}
			}		
		}
	}
	
	static void processpCorrectedEvent(DataEvent event, int corEventCounter, byte setsector, IndexedList<Double> thetaCorFac) {
		
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Track") && event.hasBank("REC::Traj"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int ecount = 0;
			int protoncount = 0;
			int e_pindex = -1;
			int proton_pindex = -1;
			float M = (float) 0.938;
			float E = (float) 6.535;
			float W = -20;
			float dtheta = -1;
			float theta_deg_e = 0;
			float phi_deg_e = -200;
			float phi_rot_deg_e = -200;
			double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0);
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0);
			Vector3D v_e = new Vector3D (0, 0, 0);
			Vector3D v_proton = new Vector3D (0, 0, 0);
			float cdchi2_proton = 100;
			short cdNDF_proton = 100;
			Vector3D dc_hit_e = new Vector3D (0, 0, 0);
			Vector3D dc_hit_rot_e = new Vector3D (0, 0, 0);
			double pcf = 0;
			for(int i = 0; i < rec.rows(); i++)
			{
				int pid = rec.getInt("pid", i);
				float px = rec.getFloat("px", i);
				float py = rec.getFloat("py", i);
				float pz = rec.getFloat("pz", i);
				float vx = rec.getFloat("vx", i);
				float vy = rec.getFloat("vy", i);
				float vz = rec.getFloat("vz", i);
				float p = (float) Math.sqrt((px*px)+(py*py)+(pz*pz));
				if(pid == 11)
				{
					ecount++;
					if(ecount == 1)
					{
						p_e = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.000511*0.000511)));
						v_e = new Vector3D (vx, vy, vz);
						e_pindex = i;
					}
				}
				if(pid == 2212)
				{
					protoncount++;
					if(protoncount == 1)
					{
						p_proton = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.938272*0.938272)));
						v_proton = new Vector3D (vx, vy, vz);
						proton_pindex = i;
					}
				}
			}
			HipoDataBank rectrac = (HipoDataBank) event.getBank("REC::Track");
			byte e_detector = 0;
			byte e_sector = 0;
			byte proton_detector = 0;
			for(int j = 0; j < rectrac.rows(); j++)
			{
				short pindex = rectrac.getShort("pindex", j);
				byte detector = rectrac.getByte("detector", j);
				byte sector = rectrac.getByte("sector", j);
				float chi2 = rectrac.getFloat("chi2", j);
				short NDF = rectrac.getShort("NDF", j);
				if(pindex == e_pindex)
				{
					e_detector = detector;
					e_sector = sector;
				}
				if(pindex == proton_pindex)
				{
					proton_detector = detector;
					cdchi2_proton = chi2;
					cdNDF_proton = NDF;
				}
			}
			HipoDataBank rectraj = (HipoDataBank) event.getBank("REC::Traj");
			for(int k = 0; k < rectraj.rows(); k++)
			{
				short pindex = rectraj.getShort("pindex", k);
				byte detector = rectraj.getByte("detector", k);
				byte layer = rectraj.getByte("layer", k);
				float x = rectraj.getFloat("x", k);
				float y = rectraj.getFloat("y", k);
				float z = rectraj.getFloat("z", k);
				if(pindex == e_pindex && detector == 6 && e_detector == 6 && layer ==  6)
				{
					dc_hit_e = new Vector3D (x, y, z);
					theta_deg_e = (float) (57.3*dc_hit_e.theta());
					phi_deg_e = (float) (57.3*dc_hit_e.phi());
					if(phi_deg_e <= -150) phi_deg_e = phi_deg_e+360;
					
					double[] b0 = new double[]{thetaCorFac.getItem(e_sector-1, 0, 0), thetaCorFac.getItem(e_sector-1, 0, 1), thetaCorFac.getItem(e_sector-1, 0, 2),
							thetaCorFac.getItem(e_sector-1, 0, 3), thetaCorFac.getItem(e_sector-1, 0, 4)};
					double[] b1 = new double[]{thetaCorFac.getItem(e_sector-1, 1, 0), thetaCorFac.getItem(e_sector-1, 1, 1), thetaCorFac.getItem(e_sector-1, 1, 2),
												thetaCorFac.getItem(e_sector-1, 1, 3), thetaCorFac.getItem(e_sector-1, 1, 4)};
					
					double a0 = (b0[4]*theta_deg_e*theta_deg_e*theta_deg_e*theta_deg_e)+(b0[3]*theta_deg_e*theta_deg_e*theta_deg_e)
									+(b0[2]*theta_deg_e*theta_deg_e)+(b0[1]*theta_deg_e)+b0[0];
					double a1 = (b1[4]*theta_deg_e*theta_deg_e*theta_deg_e*theta_deg_e)+(b1[3]*theta_deg_e*theta_deg_e*theta_deg_e)
									+(b1[2]*theta_deg_e*theta_deg_e)+(b1[1]*theta_deg_e)+b1[0];
					
					dc_hit_rot_e = new Vector3D (x*Math.cos(phi_rot[e_sector-1]/57.3)+y*Math.sin(phi_rot[e_sector-1]/57.3),
							y*Math.cos(phi_rot[e_sector-1]/57.3)-x*Math.sin(phi_rot[e_sector-1]/57.3), z);
					phi_rot_deg_e = (float) (57.3*dc_hit_rot_e.phi());
					pcf = (a1*phi_rot_deg_e)+a0;
				}
			}
			
			
			if(ecount == 1 && protoncount == 1
					&& v_e.z() > -12 && v_e.z() < 7 && Math.abs(v_e.z()-v_proton.z()) < 2.5 + (2.5/p_proton.p())
					&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30)))
				//	&& 57.3*p_proton.theta() < 75)
			{
				W = (float) Math.sqrt((M*M)+(2*M*(E-(p_e.p()*pcf)))+(2*E*p_e.p()*pcf*((p_e.pz()/p_e.p())-1)));
				dtheta = (float) ((Math.acos(1+(M/E)-(M/(p_e.p()*pcf))))/p_e.theta());
				
				hpcorWs.fill(W);
				histGroups_p_cor_W_sec.getItem(e_sector-1).fill(W);
				
				double[] theta_bnd  = new double[] {7, 8.5, 10, 12, 14, 16, 18, 21, 24, 33};
				for(int thbin = 0; thbin < 9; thbin++)
				{
					if(theta_deg_e > theta_bnd[thbin] && theta_deg_e <= theta_bnd[thbin+1])
					{
						double[] phi_rot_bnd = new double[] {-30, -9, 0, 9, 30};
						for(int phbin = 0; phbin < 4; phbin++)
						{
							if(phi_rot_deg_e > phi_rot_bnd[phbin] && phi_rot_deg_e <= phi_rot_bnd[phbin+1])
							{
								histGroups_p_cor_W_theta_phi_bin.getItem(e_sector-1, thbin, phbin).fill(W);
								histGroups_dtheta.getItem(e_sector-1, thbin, phbin).fill(dtheta);
								break;
							}
						}
						break;
					}
				}
			}
		}
	}
	
	static void processthetapCorrectedEvent(DataEvent event, int corEventCounter, byte setsector, IndexedList<Double> thetaCorFac, IndexedList<Double> thetaCorFac1) {
		
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Track") && event.hasBank("REC::Traj"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int ecount = 0;
			int protoncount = 0;
			int e_pindex = -1;
			int proton_pindex = -1;
			float M = (float) 0.938;
			float E = (float) 6.535;
			float W = -20;
			float theta_deg_e = 0;
			float phi_deg_e = -200;
			float phi_rot_deg_e = -200;
			double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0);
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0);
			Vector3D v_e = new Vector3D (0, 0, 0);
			Vector3D v_proton = new Vector3D (0, 0, 0);
			float cdchi2_proton = 100;
			short cdNDF_proton = 100;
			Vector3D dc_hit_e = new Vector3D (0, 0, 0);
			Vector3D dc_hit_rot_e = new Vector3D (0, 0, 0);
			double pcf = 0;
			double thetacf = 0;
			for(int i = 0; i < rec.rows(); i++)
			{
				int pid = rec.getInt("pid", i);
				float px = rec.getFloat("px", i);
				float py = rec.getFloat("py", i);
				float pz = rec.getFloat("pz", i);
				float vx = rec.getFloat("vx", i);
				float vy = rec.getFloat("vy", i);
				float vz = rec.getFloat("vz", i);
				float p = (float) Math.sqrt((px*px)+(py*py)+(pz*pz));
				if(pid == 11)
				{
					ecount++;
					if(ecount == 1)
					{
						p_e = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.000511*0.000511)));
						v_e = new Vector3D (vx, vy, vz);
						e_pindex = i;
					}
				}
				if(pid == 2212)
				{
					protoncount++;
					if(protoncount == 1)
					{
						p_proton = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.938272*0.938272)));
						v_proton = new Vector3D (vx, vy, vz);
						proton_pindex = i;
					}
				}
			}
			HipoDataBank rectrac = (HipoDataBank) event.getBank("REC::Track");
			byte e_detector = 0;
			byte e_sector = 0;
			byte proton_detector = 0;
			for(int j = 0; j < rectrac.rows(); j++)
			{
				short pindex = rectrac.getShort("pindex", j);
				byte detector = rectrac.getByte("detector", j);
				byte sector = rectrac.getByte("sector", j);
				float chi2 = rectrac.getFloat("chi2", j);
				short NDF = rectrac.getShort("NDF", j);
				if(pindex == e_pindex)
				{
					e_detector = detector;
					e_sector = sector;
				}
				if(pindex == proton_pindex)
				{
					proton_detector = detector;
					cdchi2_proton = chi2;
					cdNDF_proton = NDF;
				}
			}
			HipoDataBank rectraj = (HipoDataBank) event.getBank("REC::Traj");
			for(int k = 0; k < rectraj.rows(); k++)
			{
				short pindex = rectraj.getShort("pindex", k);
				byte detector = rectraj.getByte("detector", k);
				byte layer = rectraj.getByte("layer", k);
				float x = rectraj.getFloat("x", k);
				float y = rectraj.getFloat("y", k);
				float z = rectraj.getFloat("z", k);
				if(pindex == e_pindex && detector == 6 && e_detector == 6 && layer ==  6)
				{
					dc_hit_e = new Vector3D (x, y, z);
					theta_deg_e = (float) (57.3*dc_hit_e.theta());
					phi_deg_e = (float) (57.3*dc_hit_e.phi());
					if(phi_deg_e <= -150) phi_deg_e = phi_deg_e+360;
					dc_hit_rot_e = new Vector3D (x*Math.cos(phi_rot[e_sector-1]/57.3)+y*Math.sin(phi_rot[e_sector-1]/57.3),
							y*Math.cos(phi_rot[e_sector-1]/57.3)-x*Math.sin(phi_rot[e_sector-1]/57.3), z);
					phi_rot_deg_e = (float) (57.3*dc_hit_rot_e.phi());
					
					double[] b0 = new double[]{thetaCorFac.getItem(e_sector-1, 0, 0), thetaCorFac.getItem(e_sector-1, 0, 1), thetaCorFac.getItem(e_sector-1, 0, 2),
							thetaCorFac.getItem(e_sector-1, 0, 3), thetaCorFac.getItem(e_sector-1, 0, 4)};
					double[] b1 = new double[]{thetaCorFac.getItem(e_sector-1, 1, 0), thetaCorFac.getItem(e_sector-1, 1, 1), thetaCorFac.getItem(e_sector-1, 1, 2),
												thetaCorFac.getItem(e_sector-1, 1, 3), thetaCorFac.getItem(e_sector-1, 1, 4)};
					
					double a0 = (b0[4]*theta_deg_e*theta_deg_e*theta_deg_e*theta_deg_e)+(b0[3]*theta_deg_e*theta_deg_e*theta_deg_e)
									+(b0[2]*theta_deg_e*theta_deg_e)+(b0[1]*theta_deg_e)+b0[0];
					double a1 = (b1[4]*theta_deg_e*theta_deg_e*theta_deg_e*theta_deg_e)+(b1[3]*theta_deg_e*theta_deg_e*theta_deg_e)
									+(b1[2]*theta_deg_e*theta_deg_e)+(b1[1]*theta_deg_e)+b1[0];
					
					pcf = (a1*(phi_rot_deg_e))+a0;
					
					double[] b0_1 = new double[]{thetaCorFac1.getItem(e_sector-1, 0, 0), thetaCorFac1.getItem(e_sector-1, 0, 1), thetaCorFac1.getItem(e_sector-1, 0, 2),
													thetaCorFac1.getItem(e_sector-1, 0, 3), thetaCorFac1.getItem(e_sector-1, 0, 4)};
					double[] b1_1 = new double[]{thetaCorFac1.getItem(e_sector-1, 1, 0), thetaCorFac1.getItem(e_sector-1, 1, 1), thetaCorFac1.getItem(e_sector-1, 1, 2),
													thetaCorFac1.getItem(e_sector-1, 1, 3), thetaCorFac1.getItem(e_sector-1, 1, 4)};
					double[] b2_1 = new double[]{thetaCorFac1.getItem(e_sector-1, 2, 0), thetaCorFac1.getItem(e_sector-1, 2, 1), thetaCorFac1.getItem(e_sector-1, 2, 2),
													thetaCorFac1.getItem(e_sector-1, 2, 3), thetaCorFac1.getItem(e_sector-1, 2, 4)};
					
					double a0_1 = (b0_1[4]*theta_deg_e*theta_deg_e*theta_deg_e*theta_deg_e)+(b0_1[3]*theta_deg_e*theta_deg_e*theta_deg_e)+(b0_1[2]*theta_deg_e*theta_deg_e)
									+(b0_1[1]*theta_deg_e)+b0_1[0];
					double a1_1 = (b1_1[4]*theta_deg_e*theta_deg_e*theta_deg_e*theta_deg_e)+(b1_1[3]*theta_deg_e*theta_deg_e*theta_deg_e)+(b1_1[2]*theta_deg_e*theta_deg_e)
									+(b1_1[1]*theta_deg_e)+b1_1[0];
					double a2_1 = (b2_1[4]*theta_deg_e*theta_deg_e*theta_deg_e*theta_deg_e)+(b2_1[3]*theta_deg_e*theta_deg_e*theta_deg_e)+(b2_1[2]*theta_deg_e*theta_deg_e)
									+(b2_1[1]*theta_deg_e)+b2_1[0];
					
					thetacf = (a2_1*(phi_rot_deg_e)*(phi_rot_deg_e))+(a1_1*(phi_rot_deg_e))+a0_1;
				}
			}
			
			if(ecount == 1 && protoncount == 1
					&& v_e.z() > -12 && v_e.z() < 7 && Math.abs(v_e.z()-v_proton.z()) < 2.5 + (2.5/p_proton.p())
					&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30)))
				//	&& 57.3*p_proton.theta() < 75)
			{
				W = (float) Math.sqrt((M*M)+(2*M*(E-(p_e.p()*pcf)))-(4*E*p_e.p()*pcf
											*(Math.sin(p_e.theta()*thetacf/2))*(Math.sin(p_e.theta()*thetacf/2))));
				
				hthpcorWs.fill(W);
				histGroups_theta_p_cor_W_sec.getItem(e_sector-1).fill(W);
					
				double[] theta_bnd  = new double[] {7, 8.5, 10, 12, 14, 16, 18, 21, 24, 33};
				for(int thbin = 0; thbin < 9; thbin++)
				{
						if(theta_deg_e > theta_bnd[thbin] && theta_deg_e <= theta_bnd[thbin+1])
					{
							
						double[] phi_rot_bnd = new double[] {-30, -9, 0, 9, 30};
						for(int phbin = 0; phbin < 4; phbin++)
						{
							if(phi_rot_deg_e > phi_rot_bnd[phbin] && phi_rot_deg_e <= phi_rot_bnd[phbin+1])
							{
								histGroups_theta_p_cor_W_theta_phi_bin.getItem(e_sector-1, thbin, phbin).fill(W);
								break;
							}
						}
						break;
					}
				}
			}
		}
	}
	
	static class FCN implements FCNBase {
		FCN(List<Double> phibci, List<Double> meani, List<Double> sigmai)
		{
			phibc = phibci;
			mean = meani;
			sigma = sigmai;
		}
		public double errorDef()
		{
			return 1;
		}
		public double valueOf(double[] par)
		{
			double a0 = par[0];
			double a1 = par[1];
			double chisq = 0;
			for(int n = 0; n < mean.size(); n++)
			{
				double phibcn = phibc.get(n);
				double meann = mean.get(n);
				double sigman = sigma.get(n);
				double delta = meann-((a1*phibcn)+a0);
				chisq += (delta*delta/(sigman*sigman));
			}
			return chisq;
		}
		private List<Double> phibc;
		private List<Double> mean;
		private List<Double> sigma;
	}
	
	static class FCN1 implements FCNBase {
		FCN1(List<Double> thetabcfi, List<Double> cfi, List<Double> cfei)
		{
			thetabcf = thetabcfi;
			cf = cfi;
			cfe = cfei;
		}
		public double errorDef()
		{
			return 1;
		}
		public double valueOf(double[] par)
		{
			double b0 = par[0];
			double b1 = par[1];
			double b2 = par[2];
			double b3 = par[3];
			double b4 = par[4];
			double chisq = 0;
			for(int n = 0; n < cf.size(); n++)
			{
				double thetabcfn = thetabcf.get(n);
				double cfn = cf.get(n);
				double cfen = cfe.get(n);
				double delta = cfn-((b4*thetabcfn*thetabcfn*thetabcfn*thetabcfn)+(b3*thetabcfn*thetabcfn*thetabcfn)
								+(b2*thetabcfn*thetabcfn)+(b1*thetabcfn)+b0);
				chisq += (delta*delta/(cfen*cfen));
			}
			return chisq;
		}
		private List<Double> thetabcf;
		private List<Double> cf;
		private List<Double> cfe;
	}
	
	static class FCN2 implements FCNBase {
		FCN2(List<Double> phibci, List<Double> mean1i, List<Double> sigma1i)
		{
			phibc = phibci;
			mean1 = mean1i;
			sigma1 = sigma1i;
		}
		public double errorDef()
		{
			return 1;
		}
		public double valueOf(double[] par)
		{
			double a0 = par[0];
			double a1 = par[1];
			double a2 = par[2];
			double chisq = 0;
			for(int n = 0; n < mean1.size(); n++)
			{
				double phibcn = phibc.get(n);
				double mean1n = mean1.get(n);
				double sigma1n = sigma1.get(n);
				double delta = mean1n-((a2*phibcn*phibcn)+(a1*phibcn)+a0);;
				chisq += (delta*delta/(sigma1n*sigma1n));
			}
			return chisq;
		}
		private List<Double> phibc;
		private List<Double> mean1;
		private List<Double> sigma1;
	}
	
	static class FCN3 implements FCNBase {
		FCN3(List<Double> thetabcfi, List<Double> cf1i, List<Double> cfe1i)
		{
			thetabcf = thetabcfi;
			cf1 = cf1i;
			cfe1 = cfe1i;
		}
		public double errorDef()
		{
			return 1;
		}
		public double valueOf(double[] par)
		{
			double b0 = par[0];
			double b1 = par[1];
			double b2 = par[2];
			double b3 = par[3];
			double b4 = par[4];
			double chisq = 0;
			for(int n = 0; n < cf1.size(); n++)
			{
				double thetabcfn = thetabcf.get(n);
				double cf1n = cf1.get(n);
				double cfe1n = cfe1.get(n);
				double delta = cf1n-((b4*thetabcfn*thetabcfn*thetabcfn*thetabcfn)+(b3*thetabcfn*thetabcfn*thetabcfn)
								+(b2*thetabcfn*thetabcfn)+(b1*thetabcfn)+b0);
				chisq += (delta*delta/(cfe1n*cfe1n));
			}
			return chisq;
		}
		private List<Double> thetabcf;
		private List<Double> cf1;
		private List<Double> cfe1;
	}
	
	public static void main(String[] args) {
		
		histos_W_sec();
		histos_p_cor_W_sec();
		histos_theta_p_cor_W_sec();
		histos2();
		histos_W_vs_phi_bin();
		histos_p_cor_W_vs_phi_bin();
		histos_theta_p_cor_W_theta_phi_bin();
		histos_dp();
		histos_dtheta();
		
		byte setsector = 1;
			
			reader.open("C:/Users/joshtanj/Documents/download/merged_skim_e_bank_6535MeV_skim_elastic_5950_5899_5963_5958_5885.hipo");
			
			
			int eventCounter = 0;
			while(reader.hasEvent()/* && eventCounter < 5000000*/)
			{
				eventCounter++;
				processEvent(reader.getNextEvent(), eventCounter, setsector);
				if(eventCounter%50000 == 0) System.out.println("Event: " + eventCounter);
			}
			
			reader.close();
			
			IndexedList<Double> Wmean = new IndexedList<>(2);
			IndexedList<Double> Wsigma = new IndexedList<>(2);
			
			JFrame frameW = new JFrame("W");
			frameW.setSize(1500, 1000);
			EmbeddedCanvas canhg1 = new EmbeddedCanvas();
			frameW.add(canhg1);
			frameW.setLocationRelativeTo(null);
			frameW.setVisible(true);
			canhg1.divide(3, 2);
			for(int canhg1i = 0; canhg1i < 6; canhg1i++){
				canhg1.cd(canhg1i);
				canhg1.setFont("Arial");
				histGroups_W_sec.getItem(canhg1i).setTitle("Sector " + (canhg1i+1) + " W");
				histGroups_W_sec.getItem(canhg1i).setTitleX("W [GeV]");
				histGroups_W_sec.getItem(canhg1i).setTitleY("Counts");
				histGroups_W_sec.getItem(canhg1i).setOptStat(10);
				canhg1.getPad(canhg1i).setTitleFontSize(20);
				canhg1.getPad(canhg1i).setAxisTitleFontSize(20);
				canhg1.getPad(canhg1i).setAxisLabelFontSize(20);
				canhg1.getPad(canhg1i).setStatBoxFontSize(20);
				canhg1.draw(histGroups_W_sec.getItem(canhg1i), "same");
				F1D f001 = new F1D("f001", "[amp]*gaus(x,[mean],[sigma])",
									0.7+(histGroups_W_sec.getItem(canhg1i).getMaximumBin()*(1.2-0.7)/200)-0.035,
									0.7+(histGroups_W_sec.getItem(canhg1i).getMaximumBin()*(1.2-0.7)/200)+0.025);
				f001.setParameter(0, histGroups_W_sec.getItem(canhg1i).getMax()*1.2);
				f001.setParameter(1, 0.7+(histGroups_W_sec.getItem(canhg1i).getMaximumBin()*(1.2-0.7)/200));
				f001.setParameter(2, 0.01);
				DataFitter.fit(f001, histGroups_W_sec.getItem(canhg1i), "Q");
				f001.setLineColor(2);
				f001.setLineWidth(3);
				f001.setOptStat(11110);
				canhg1.draw(f001, "same");
				Wmean.add(f001.getParameter(1), 0, canhg1i);
				Wsigma.add(f001.getParameter(2), 0, canhg1i);
			}
			
			JFrame frameWsec1 = new JFrame("All Sector W");
			frameWsec1.setSize(500, 500);
			EmbeddedCanvas canhs1 = new EmbeddedCanvas();
			frameWsec1.add(canhs1);
			frameWsec1.setLocationRelativeTo(null);
			frameWsec1.setVisible(true);
			canhs1.divide(1, 1);
			canhs1.cd(0);
			canhs1.setFont("Arial");
			hWs.setTitle("All Sector W");
			hWs.setTitleX("W [GeV]");
			hWs.setTitleY("Counts");
			hWs.setOptStat(10);
			canhs1.getPad(0).setTitleFontSize(20);
			canhs1.getPad(0).setAxisTitleFontSize(20);
			canhs1.getPad(0).setAxisLabelFontSize(20);
			canhs1.getPad(0).setStatBoxFontSize(20);
			canhs1.draw(hWs, "same");
			F1D f00s1 = new F1D("f00s1", "[amp]*gaus(x,[mean],[sigma])",
								0.7+(hWs.getMaximumBin()*(1.2-0.7)/200)-0.035,
								0.7+(hWs.getMaximumBin()*(1.2-0.7)/200)+0.025);
			f00s1.setParameter(0, hWs.getMax()*1.2);
			f00s1.setParameter(1, 0.7+(hWs.getMaximumBin()*(1.2-0.7)/200));
			f00s1.setParameter(2, 0.01);
			DataFitter.fit(f00s1, hWs, "Q");
			f00s1.setLineColor(2);
			f00s1.setLineWidth(3);
			f00s1.setOptStat(11110);
			canhs1.draw(f00s1, "same");

			double[] theta_bnd  = new double[] {7, 8.5, 10, 12, 14, 16, 18, 21, 24, 33};
			
			for(int canhg2i = 0; canhg2i < 9; canhg2i++){
				JFrame frameWphi = new JFrame("W vs. #phi");
				frameWphi.setSize(1500, 500);
				EmbeddedCanvas canhg2 = new EmbeddedCanvas();
				frameWphi.add(canhg2);
				frameWphi.setLocationRelativeTo(null);
				frameWphi.setVisible(true);
				canhg2.divide(1, 1);
				canhg2.cd(0);
				canhg2.setFont("Arial");
				histGroups_W_vs_phi.getItem(canhg2i).setTitle(theta_bnd[canhg2i] + " #degree < #theta <= "
																+ theta_bnd[canhg2i+1] + " #degree W vs. #phi");
				histGroups_W_vs_phi.getItem(canhg2i).setTitleX("#phi [#degree]");
				histGroups_W_vs_phi.getItem(canhg2i).setTitleY("W [GeV]");
				canhg2.getPad(0).setTitleFontSize(20);
				canhg2.getPad(0).setAxisTitleFontSize(20);
				canhg2.getPad(0).setAxisLabelFontSize(20);
				canhg2.draw(histGroups_W_vs_phi.getItem(canhg2i));
			}
			

			double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
			double[] phi_rot_bnd = new double[] {-30, -9, 0, 9, 30};
			
			IndexedList<Double> mean = new IndexedList<>(3);
			IndexedList<Double> sigma = new IndexedList<>(3);
			
			IndexedList<Double> meanncor = new IndexedList<>(3);
			IndexedList<Double> sigmancor = new IndexedList<>(3);
			
			for(int seci = 0; seci < 6; seci++)
			{
				for(int theta_bini = 0; theta_bini < 9; theta_bini++)
				{
					for(int phi_bini = 0; phi_bini < 4; phi_bini++)
					{
						F1D f003 = new F1D("f003", "[amp]*gaus(x,[mean],[sigma])", 
											0.7+(histGroups_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*(1.2-0.7)/200)-0.065,
											0.7+(histGroups_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*(1.2-0.7)/200)+0.045);
						f003.setParameter(0, histGroups_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax());
						f003.setParameter(1, 0.7+(histGroups_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*(1.2-0.7)/200));
						f003.setParameter(2, 0.01);
						DataFitter.fit(f003, histGroups_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini), "Q");
						meanncor.add(f003.getParameter(1), seci, theta_bini, phi_bini);
						sigmancor.add(f003.getParameter(2), seci, theta_bini, phi_bini);
					}
				}
			}
			
			for(int theta_bini = 0; theta_bini < 9; theta_bini++)
			{
				JFrame framedpbin = new JFrame("#deltap Binned");
				framedpbin.setSize(1000, 1000);
				EmbeddedCanvas canhgdp = new EmbeddedCanvas();
				framedpbin.add(canhgdp);
				framedpbin.setLocationRelativeTo(null);
				framedpbin.setVisible(true);
				canhgdp.divide(2, 2);
				for(int phi_bini = 0; phi_bini < 4; phi_bini++)
				{
					canhgdp.cd(phi_bini);
					canhgdp.setFont("Arial");
					histGroups_dp.getItem(setsector-1, theta_bini, phi_bini).setTitle("Sector " + setsector + " e^-: " + theta_bnd[theta_bini]
																						+ " #degree < #theta <= " + theta_bnd[theta_bini+1] + " #degree, "
																						+ (phi_rot[setsector-1]+phi_rot_bnd[phi_bini])
																						+ " #degree < #phi <= "
																						+ (phi_rot[setsector-1]+phi_rot_bnd[phi_bini+1]) + " #degree");
					histGroups_dp.getItem(setsector-1, theta_bini, phi_bini).setTitleX("#deltap");
					histGroups_dp.getItem(setsector-1, theta_bini, phi_bini).setTitleY("Counts");
					histGroups_dp.getItem(setsector-1, theta_bini, phi_bini).setOptStat(10);
					canhgdp.getPad(phi_bini).setTitleFontSize(20);
					canhgdp.getPad(phi_bini).setAxisTitleFontSize(20);
					canhgdp.getPad(phi_bini).setAxisLabelFontSize(20);
					canhgdp.getPad(phi_bini).setStatBoxFontSize(20);
					canhgdp.draw(histGroups_dp.getItem(setsector-1, theta_bini, phi_bini));
					F1D f00dp = new F1D("f00dp", "[amp]*gaus(x,[mean],[sigma])",
										0.85+(histGroups_dp.getItem(setsector-1, theta_bini, phi_bini).getMaximumBin()*(1.15-0.85)/200)-0.0175,
										0.85+(histGroups_dp.getItem(setsector-1, theta_bini, phi_bini).getMaximumBin()*(1.15-0.85)/200)+0.0125);
					f00dp.setParameter(0, histGroups_dp.getItem(setsector-1, theta_bini, phi_bini).getMax());
					f00dp.setParameter(1, 0.85+(histGroups_dp.getItem(setsector-1, theta_bini, phi_bini).getMaximumBin()*(1.15-0.85)/200));
					f00dp.setParameter(2, 0.01);
					DataFitter.fit(f00dp, histGroups_dp.getItem(setsector-1, theta_bini, phi_bini), "Q");
					f00dp.setLineColor(2);
					f00dp.setLineWidth(3);
					f00dp.setOptStat(11110);
					canhgdp.draw(f00dp, "same");
				}
			}
			
			for(int secdpi = 0; secdpi < 6; secdpi++)
			{
				for(int theta_bini = 0; theta_bini < 9; theta_bini++)
				{
					for(int phi_bini = 0; phi_bini < 4; phi_bini++)
					{
						F1D f00dp = new F1D("f00dp", "[amp]*gaus(x,[mean],[sigma])",
											0.85+(histGroups_dp.getItem(secdpi, theta_bini, phi_bini).getMaximumBin()*(1.15-0.85)/200)-0.0175,
											0.85+(histGroups_dp.getItem(secdpi, theta_bini, phi_bini).getMaximumBin()*(1.15-0.85)/200)+0.0125);
						f00dp.setParameter(0, histGroups_dp.getItem(secdpi, theta_bini, phi_bini).getMax());
						f00dp.setParameter(1, 0.85+(histGroups_dp.getItem(secdpi, theta_bini, phi_bini).getMaximumBin()*(1.15-0.85)/200));
						f00dp.setParameter(2, 0.01);
						DataFitter.fit(f00dp, histGroups_dp.getItem(secdpi, theta_bini, phi_bini), "Q");
						mean.add(f00dp.getParameter(1), secdpi, theta_bini, phi_bini);
						sigma.add(f00dp.parameter(1).error(), secdpi, theta_bini, phi_bini);
					}
				}
			}
			
			double[] phirotba = new double[] {-13.27, -4.05, 4.05, 13.27};
			
			IndexedList<GraphErrors> histerrGroups_W_vs_phi = new IndexedList<GraphErrors>(2);
			
			for(int ihisterrsector = 0; ihisterrsector < 6; ihisterrsector++) {
				for(int itheta_bin = 0; itheta_bin < 9; itheta_bin++) {
					GraphErrors heW_vs_phi = new GraphErrors();
					histerrGroups_W_vs_phi.add(heW_vs_phi, ihisterrsector, itheta_bin);
				}
			}
			
			for(int seci = 0; seci < 6; seci++)
			{
				for(int thetai = 0; thetai < 9; thetai++)
				{
					for(int phii = 0; phii < 4; phii++)
					{
						histerrGroups_W_vs_phi.getItem(seci, thetai).addPoint((phirotba[phii]+phi_rot[seci]), meanncor.getItem(seci, thetai,  phii),
																				0, sigmancor.getItem(seci, thetai, phii));
					}
				}
			}
			
			for(int thetai = 0; thetai < 9; thetai++){
				JFrame frameWphibin = new JFrame("W vs. #phi Binned");
				frameWphibin.setSize(1500, 1000);
				EmbeddedCanvas canhge = new EmbeddedCanvas();
				frameWphibin.add(canhge);
				frameWphibin.setLocationRelativeTo(null);
				frameWphibin.setVisible(true);
				canhge.divide(3, 2);
				for(int seci = 0; seci < 6; seci++){
					canhge.cd(seci);
					canhge.setFont("Arial");
					histerrGroups_W_vs_phi.getItem(seci, thetai).setTitle("Sector " + (seci+1) + " e^-: " + theta_bnd[thetai] + " #degree < #theta <= "
																			+ theta_bnd[thetai+1] + " #degree, W vs. #phi");
					histerrGroups_W_vs_phi.getItem(seci, thetai).setTitleX("#phi [#degree]");
					histerrGroups_W_vs_phi.getItem(seci, thetai).setTitleY("W [GeV]");
					histerrGroups_W_vs_phi.getItem(seci, thetai).setMarkerSize(5);
					histerrGroups_W_vs_phi.getItem(seci, thetai).setMarkerStyle(1);
					histerrGroups_W_vs_phi.getItem(seci, thetai).setLineThickness(1);
					canhge.getPad(seci).setTitleFontSize(20);
					canhge.getPad(seci).setAxisTitleFontSize(20);
					canhge.getPad(seci).setAxisLabelFontSize(20);
					canhge.getPad(seci).setAxisRange((phi_rot[seci]-30), (phi_rot[seci]+30), 0.85, 1.15);
					canhge.draw(histerrGroups_W_vs_phi.getItem(seci, thetai));
					}
			}
			
			double[] thetamean = new double[]{7.62, 9.12, 10.79, 12.80, 14.83, 16.82, 19.47, 22.28, 26.24};
			
			IndexedList<Double> phiCorFac = new IndexedList<>(3);
			
			for(int seci = 0; seci < 6; seci++)
			{
				for(int thetabi = 0; thetabi < 9; thetabi++)
				{
					List<Double> meani = new ArrayList<>();
					List<Double> sigmai = new ArrayList<>();
					List<Double> phibci = new ArrayList<>();
					
					for(int phibi = 0; phibi < 4; phibi++)
					{
						if(mean.getItem(seci, thetabi, phibi) > 0.95 && mean.getItem(seci, thetabi, phibi) < 1.025
								&& Math.abs(sigma.getItem(seci, thetabi, phibi)) < 0.01)
						{
							meani.add(mean.getItem(seci, thetabi, phibi));
							sigmai.add(sigma.getItem(seci, thetabi, phibi));
							phibci.add(phirotba[phibi]);
						}
					}
					
					FCN theFCN = new FCN(phibci, meani, sigmai);
						
					MnUserParameters upar = new MnUserParameters();
					upar.add("a0", 0.94, 0.1);
					upar.add("a1", 0, 0.01);
					System.out.println("Initial parameters: "+upar); 
				    
				    System.out.println("start migrad");
					MnMigrad migrad = new MnMigrad(theFCN, upar); 
				    FunctionMinimum min = migrad.minimize(); 
				    if(!min.isValid()) 
				    { 
				    	//try with higher strategy 
					    	System.out.println("FM is invalid, try with strategy = 2."); 
					    	MnMigrad migrad2 = new MnMigrad(theFCN, min.userState(), new MnStrategy(2)); 
					    	min = migrad2.minimize(); 
				    }
				    System.out.println("minimum: "+min);
				    for(int cfi = 0; cfi < 2; cfi++)
				    {
				    	phiCorFac.add(min.userParameters().value(cfi), seci, thetabi, cfi);
				    }
				    for(int cfei = 2; cfei < 4; cfei++)
				    {
				    	phiCorFac.add(min.userParameters().error(cfei-2), seci, thetabi, cfei);
				    }
				    phiCorFac.add(Math.sqrt(min.fval()), seci, thetabi, 4);
				}
			}
			
			IndexedList<GraphErrors> histerrGroups_dp_vs_phi = new IndexedList<GraphErrors>(2);
			
			for(int ihisterrsectordp = 0; ihisterrsectordp < 6; ihisterrsectordp++) {
				for(int itheta_bindp = 0; itheta_bindp < 9; itheta_bindp++) {
					GraphErrors hedp_vs_phi = new GraphErrors();
					histerrGroups_dp_vs_phi.add(hedp_vs_phi, ihisterrsectordp, itheta_bindp);
				}
			}
			
			for(int seci = 0; seci < 6; seci++)
			{
				for(int thetadpi = 0; thetadpi < 9; thetadpi++)
				{
					for(int phidpi = 0; phidpi < 4; phidpi++)
					{
						histerrGroups_dp_vs_phi.getItem(seci, thetadpi).addPoint(phirotba[phidpi], mean.getItem(seci, thetadpi,  phidpi),
																					0, sigma.getItem(seci, thetadpi, phidpi));
					}
				}
			}
			
			for(int framehgedpi = 0; framehgedpi < 9; framehgedpi++){
				JFrame framedpphibin = new JFrame("#deltap vs. #phi Binned");
				framedpphibin.setSize(1500, 1000);
				EmbeddedCanvas canhgedp = new EmbeddedCanvas();
				framedpphibin.add(canhgedp);
				framedpphibin.setLocationRelativeTo(null);
				framedpphibin.setVisible(true);
				canhgedp.divide(3, 2);
				for(int seci = 0; seci < 6; seci++)
				{
					canhgedp.cd(seci);
					canhgedp.setFont("Arial");
					histerrGroups_dp_vs_phi.getItem(seci, framehgedpi).setTitle("Sector " + (seci+1) + " e^-: "
																					+ theta_bnd[framehgedpi]
																					+ " #degree < #theta <= "
																					+ theta_bnd[framehgedpi+1] + " #degree, "
																					+ " #deltap vs. #phi");
					histerrGroups_dp_vs_phi.getItem(seci, framehgedpi).setTitleX("#phi_midsector [#degree]");
					histerrGroups_dp_vs_phi.getItem(seci, framehgedpi).setTitleY("#deltap");
					histerrGroups_dp_vs_phi.getItem(seci, framehgedpi).setMarkerSize(5);
					histerrGroups_dp_vs_phi.getItem(seci, framehgedpi).setMarkerStyle(1);
					histerrGroups_dp_vs_phi.getItem(seci, framehgedpi).setLineThickness(1);
					canhgedp.getPad(seci).setTitleFontSize(20);
					canhgedp.getPad(seci).setAxisTitleFontSize(20);
					canhgedp.getPad(seci).setAxisLabelFontSize(20);
					canhgedp.getPad(seci).setAxisRange(-30, 30, 0.95, 1.025);
					canhgedp.draw(histerrGroups_dp_vs_phi.getItem(seci, framehgedpi));//, "L");
					F1D fdp = new F1D("fdp", "[a0]+([a1]*x)", -30, +30);
					fdp.setParameter(0, phiCorFac.getItem(seci, framehgedpi, 0));
					fdp.setParameter(1, phiCorFac.getItem(seci, framehgedpi, 1));
					fdp.setLineColor(2);
					fdp.setLineWidth(3);
					fdp.setOptStat(1110);
					canhgedp.draw(fdp, "same");
				}
			}
			
			IndexedList<Double> thetaCorFac = new IndexedList<>(3);
			
			for(int seci = 0; seci < 6; seci++)
			{
				for(int phcfi = 0; phcfi < 2; phcfi++)
				{
					List<Double> thetabcfi = new ArrayList<>();
					List<Double> cfi = new ArrayList<>();
					List<Double> cfei = new ArrayList<>();
					
					for(int thetabfi = 0; thetabfi < 9; thetabfi++)
					{
						{
							thetabcfi.add(thetamean[thetabfi]);
							cfi.add(phiCorFac.getItem(seci, thetabfi, phcfi));
							cfei.add(phiCorFac.getItem(seci, thetabfi, phcfi+2));
						}
					}
					
					FCN1 theFCN1 = new FCN1(thetabcfi, cfi, cfei);
						
					MnUserParameters upar1 = new MnUserParameters();
					upar1.add("b0", 0, 0.1);
					upar1.add("b1", 0, 0.01);
					upar1.add("b2", 0, 0.001);
					upar1.add("b3", 0, 0.0001);
					upar1.add("b4", 0, 0.00001);
					System.out.println("Initial parameters: "+upar1); 
				    
				    System.out.println("start migrad");
					MnMigrad migrad1 = new MnMigrad(theFCN1, upar1); 
				    FunctionMinimum min1 = migrad1.minimize(); 
				    if(!min1.isValid()) 
				    { 
				    	//try with higher strategy 
					    	System.out.println("FM is invalid, try with strategy = 2."); 
					    	MnMigrad migrad21 = new MnMigrad(theFCN1, min1.userState(), new MnStrategy(2)); 
					    	min1 = migrad21.minimize(); 
				    }
				    System.out.println("minimum: "+min1);
				    
				    for(int thcfi = 0; thcfi < 5; thcfi++)
				    {
				    	thetaCorFac.add(min1.userParameters().value(thcfi), seci, phcfi, thcfi);
				    }
				    for(int thcfei = 5; thcfei < 10; thcfei++)
				    {
				    	thetaCorFac.add(min1.userParameters().error(thcfei-5), seci, phcfi, thcfei);
				    }
				    thetaCorFac.add(Math.sqrt(min1.fval()), seci, phcfi, 10);
				}
			}
			
			IndexedList<GraphErrors> histerrGroupsCorFac_p = new IndexedList<GraphErrors>(2);
			
			for(int ihisterrsectorcorfac = 0; ihisterrsectorcorfac < 6; ihisterrsectorcorfac++) {
				for(int itheta_bincorfac = 0; itheta_bincorfac < 5; itheta_bincorfac++) {
					GraphErrors heWcorfac = new GraphErrors();
					histerrGroupsCorFac_p.add(heWcorfac, ihisterrsectorcorfac, itheta_bincorfac);
				}
			}			
			
			for(int seci = 0; seci < 6; seci++)
			{
				for(int cci = 0; cci < 2; cci++)
				{
					GraphErrors hhecf = histerrGroupsCorFac_p.getItem(seci, cci);
					for(int thetabgci = 0; thetabgci < 9; thetabgci++)
					{
						hhecf.addPoint(thetamean[thetabgci], phiCorFac.getItem(seci, thetabgci, cci),
										0, phiCorFac.getItem(seci, thetabgci, cci+2));
					}
				}
			}
			
			for(int framehgecfi = 0; framehgecfi < 2; framehgecfi++){
				JFrame frameWphibincorfac = new JFrame("#deltap: #phi Correction Coefficients vs. #theta Binned");
				frameWphibincorfac.setSize(1500, 1000);
				EmbeddedCanvas canhgecf = new EmbeddedCanvas();
				frameWphibincorfac.add(canhgecf);
				frameWphibincorfac.setLocationRelativeTo(null);
				frameWphibincorfac.setVisible(true);
				canhgecf.divide(3, 2);
				for(int seci = 0; seci < 6; seci++)
				{
					canhgecf.cd(seci);
					canhgecf.setFont("Arial");
					histerrGroupsCorFac_p.getItem(seci, framehgecfi).setTitle("Sector " + (seci+1) + " e^- #deltap: #phi Correction Coefficient a"
																				+ framehgecfi + " vs. #theta");
					histerrGroupsCorFac_p.getItem(seci, framehgecfi).setTitleX("#theta [#degree]");
					histerrGroupsCorFac_p.getItem(seci, framehgecfi).setTitleY("a" + framehgecfi + " [arb. unit]");
					histerrGroupsCorFac_p.getItem(seci, framehgecfi).setMarkerSize(5);
					histerrGroupsCorFac_p.getItem(seci, framehgecfi).setMarkerStyle(1);
					histerrGroupsCorFac_p.getItem(seci, framehgecfi).setLineThickness(1);
					canhgecf.getPad(seci).setTitleFontSize(20);
					canhgecf.getPad(seci).setAxisTitleFontSize(20);
					canhgecf.getPad(seci).setAxisLabelFontSize(20);
					if(framehgecfi == 0) canhgecf.getPad(seci).setAxisRange(7, 33, 0.995, 1.01);
					if(framehgecfi == 1) canhgecf.getPad(seci).setAxisRange(7, 33, -0.0005, 0.0007);
					canhgecf.draw(histerrGroupsCorFac_p.getItem(seci, framehgecfi), "same");
					F1D fcf = new F1D("fcf", "[b0]+([b1]*x)+([b2]*x*x)+([b3]*x*x*x)+([b4]*x*x*x*x)", 7, 33);
					fcf.setParameter(0, thetaCorFac.getItem(seci, framehgecfi, 0));
					fcf.setParameter(1, thetaCorFac.getItem(seci, framehgecfi, 1));
					fcf.setParameter(2, thetaCorFac.getItem(seci, framehgecfi, 2));
					fcf.setParameter(3, thetaCorFac.getItem(seci, framehgecfi, 3));
					fcf.setParameter(4, thetaCorFac.getItem(seci, framehgecfi, 4));
					fcf.setLineColor(2);
					fcf.setLineWidth(3);
					fcf.setOptStat(111110);
					canhgecf.draw(fcf, "same");
				}
			}
			
			readerc.open("C:/Users/joshtanj/Documents/download/merged_skim_e_bank_6535MeV_skim_elastic_5950_5899_5963_5958_5885.hipo");
			
			int corEventCounter = 0;
			while(readerc.hasEvent())// && corEventCounter < 500000)
			{
				corEventCounter++;
				processpCorrectedEvent(readerc.getNextEvent(), corEventCounter, setsector, thetaCorFac);
				if(corEventCounter%50000 == 0) System.out.println("Corrected Event: " + corEventCounter);
			}
			
			readerc.close();
			
			JFrame frameWcor = new JFrame("p Corrected W");
			frameWcor.setSize(1500, 1000);
			EmbeddedCanvas canhg1cor = new EmbeddedCanvas();
			frameWcor.add(canhg1cor);
			frameWcor.setLocationRelativeTo(null);
			frameWcor.setVisible(true);
			canhg1cor.divide(3, 2);
			for(int canhg1cori = 0; canhg1cori < 6; canhg1cori++){
				canhg1cor.cd(canhg1cori);
				canhg1cor.setFont("Arial");
				histGroups_p_cor_W_sec.getItem(canhg1cori).setTitle("Sector " + (canhg1cori+1) + " W");
				histGroups_p_cor_W_sec.getItem(canhg1cori).setTitleX("W [GeV]");
				histGroups_p_cor_W_sec.getItem(canhg1cori).setTitleY("Counts");
				histGroups_p_cor_W_sec.getItem(canhg1cori).setOptStat(10);
				canhg1cor.getPad(canhg1cori).setTitleFontSize(20);
				canhg1cor.getPad(canhg1cori).setAxisTitleFontSize(20);
				canhg1cor.getPad(canhg1cori).setAxisLabelFontSize(20);
				canhg1cor.getPad(canhg1cori).setStatBoxFontSize(20);
				canhg1cor.draw(histGroups_p_cor_W_sec.getItem(canhg1cori), "same");
				F1D f001c = new F1D("f001c", "[amp]*gaus(x,[mean],[sigma])",
									0.7+(histGroups_p_cor_W_sec.getItem(canhg1cori).getMaximumBin()*(1.2-0.7)/200)-0.030,
									0.7+(histGroups_p_cor_W_sec.getItem(canhg1cori).getMaximumBin()*(1.2-0.7)/200)+0.020);
				f001c.setParameter(0, histGroups_p_cor_W_sec.getItem(canhg1cori).getMax()*1.2);
				f001c.setParameter(1, 0.7+(histGroups_p_cor_W_sec.getItem(canhg1cori).getMaximumBin()*(1.2-0.7)/200));
				f001c.setParameter(2, 0.01);
				DataFitter.fit(f001c, histGroups_p_cor_W_sec.getItem(canhg1cori), "Q");
				f001c.setLineColor(2);
				f001c.setLineWidth(3);
				f001c.setOptStat(11110);
				canhg1cor.draw(f001c, "same");
				Wmean.add(f001c.getParameter(1), 1, canhg1cori);
				Wsigma.add(f001c.getParameter(2), 1, canhg1cori);
			}
			
			JFrame framecorWsec1 = new JFrame("All Sector p Corrected W");
			framecorWsec1.setSize(500, 500);
			EmbeddedCanvas canhs1c = new EmbeddedCanvas();
			framecorWsec1.add(canhs1c);
			framecorWsec1.setLocationRelativeTo(null);
			framecorWsec1.setVisible(true);
			canhs1c.divide(1, 1);
			canhs1c.cd(0);
			canhs1c.setFont("Arial");
			hpcorWs.setTitle("All Sector p Corrected W");
			hpcorWs.setTitleX("W [GeV]");
			hpcorWs.setTitleY("Counts");
			hpcorWs.setOptStat(10);
			canhs1c.getPad(0).setTitleFontSize(20);
			canhs1c.getPad(0).setAxisTitleFontSize(20);
			canhs1c.getPad(0).setAxisLabelFontSize(20);
			canhs1c.getPad(0).setStatBoxFontSize(20);
			canhs1c.draw(hpcorWs, "same");
			F1D f001c = new F1D("f001c", "[amp]*gaus(x,[mean],[sigma])",
								0.7+(hpcorWs.getMaximumBin()*(1.2-0.7)/200)-0.030,
								0.7+(hpcorWs.getMaximumBin()*(1.2-0.7)/200)+0.020);
			f001c.setParameter(0, hpcorWs.getMax()*1.2);
			f001c.setParameter(1, 0.7+(hpcorWs.getMaximumBin()*(1.2-0.7)/200));
			f001c.setParameter(2, 0.01);
			DataFitter.fit(f001c, hpcorWs, "Q");
			f001c.setLineColor(2);
			f001c.setLineWidth(3);
			f001c.setOptStat(11110);
			canhs1c.draw(f001c, "same");
			
			IndexedList<Double> meancor = new IndexedList<>(3);
			IndexedList<Double> sigmacor = new IndexedList<>(3);
			
			for(int seci = 0; seci < 6; seci++)
			{
				for(int theta_binci = 0; theta_binci < 9; theta_binci++)
				{
					for(int phi_binci = 0; phi_binci < 4; phi_binci++)
					{
						F1D f003c = new F1D("f003c", "[amp]*gaus(x,[mean],[sigma])",
											0.7+(histGroups_p_cor_W_theta_phi_bin.getItem(seci, theta_binci, phi_binci).getMaximumBin()*(1.2-0.7)/200)-0.055,
											0.7+(histGroups_p_cor_W_theta_phi_bin.getItem(seci, theta_binci, phi_binci).getMaximumBin()*(1.2-0.7)/200)+0.035);
						f003c.setParameter(0, histGroups_p_cor_W_theta_phi_bin.getItem(seci, theta_binci, phi_binci).getMax()*1.2);
						f003c.setParameter(1, 0.7+(histGroups_p_cor_W_theta_phi_bin.getItem(seci, theta_binci, phi_binci).getMaximumBin()*(1.2-0.7)/200));
						f003c.setParameter(2, 0.01);
						DataFitter.fit(f003c, histGroups_p_cor_W_theta_phi_bin.getItem(seci, theta_binci, phi_binci), "Q");
						System.out.println(f003c.getParameter(1) + "   " + f003c.getParameter(2));
						meancor.add(f003c.getParameter(1), seci, theta_binci, phi_binci);
						sigmacor.add(f003c.getParameter(2), seci,  theta_binci, phi_binci);
					}
				}
			}
			
			IndexedList<Double> mean1 = new IndexedList<>(3);
			IndexedList<Double> sigma1 = new IndexedList<>(3);
			
			for(int secdthi = 0; secdthi < 6; secdthi++)
			{
				for(int theta_bini = 0; theta_bini < 9; theta_bini++)
				{
					for(int phi_bini = 0; phi_bini < 4; phi_bini++)
					{
						F1D f00dth = new F1D("f00dth", "[amp]*gaus(x,[mean],[sigma])",
												(histGroups_dtheta.getItem(secdthi, theta_bini, phi_bini).getMaximumBin()*2.00/200)-0.15,
												(histGroups_dtheta.getItem(secdthi, theta_bini, phi_bini).getMaximumBin()*2.00/200)+0.15);
						f00dth.setParameter(0, histGroups_dtheta.getItem(secdthi, theta_bini, phi_bini).getMax());
						f00dth.setParameter(1, histGroups_dtheta.getItem(secdthi, theta_bini, phi_bini).getMaximumBin()*2.00/200);
						f00dth.setParameter(2, 0.07);
						DataFitter.fit(f00dth, histGroups_dtheta.getItem(secdthi, theta_bini, phi_bini), "Q");
						System.out.println(f00dth.getParameter(1) + "   " + f00dth.getParameter(2));
						mean1.add(f00dth.getParameter(1), secdthi, theta_bini, phi_bini);
						sigma1.add(f00dth.parameter(1).error(), secdthi, theta_bini, phi_bini);
					}
				}
			}
			
			IndexedList<Double> phiCorFac1 = new IndexedList<>(3);
			
			for(int seci = 0; seci < 6; seci++)
			{
				for(int thetabi = 0; thetabi < 9; thetabi++)
				{
					List<Double> mean1i = new ArrayList<>();
					List<Double> sigma1i = new ArrayList<>();
					List<Double> phibci = new ArrayList<>();
					
					for(int phibi = 0; phibi < 4; phibi++)
					{
						if(mean1.getItem(seci, thetabi, phibi) > 0.85 && mean1.getItem(seci, thetabi, phibi) < 1.15 && sigma1.getItem(seci, thetabi, phibi) < 0.1)
						{
							mean1i.add(mean1.getItem(seci, thetabi, phibi));
							sigma1i.add(sigma1.getItem(seci, thetabi, phibi));
							phibci.add(phi_rot[seci]);
						}
					}
					
					FCN2 theFCN2 = new FCN2(phibci, mean1i, sigma1i);
						
					MnUserParameters upar = new MnUserParameters();
					upar.add("a0", 1, 0.1);
					upar.add("a1", 0, 0.01);
					upar.add("a2", 0, 0.001);
					System.out.println("Initial parameters: "+upar); 
				    
				    System.out.println("start migrad");
					MnMigrad migrad = new MnMigrad(theFCN2, upar); 
				    FunctionMinimum min = migrad.minimize();
				    if(!min.isValid()) 
				    { 
				    	//try with higher strategy 
					    	System.out.println("FM is invalid, try with strategy = 2."); 
					    	MnMigrad migrad2 = new MnMigrad(theFCN2, min.userState(), new MnStrategy(2)); 
					    	min = migrad2.minimize(); 
				    }
				    System.out.println("minimum: "+min);
				    for(int cfi = 0; cfi < 3; cfi++)
				    {
				    	phiCorFac1.add(min.userParameters().value(cfi), seci, thetabi, cfi);
				    }
				    for(int cfei = 3; cfei < 6; cfei++)
				    {
				    	phiCorFac1.add(min.userParameters().error(cfei-3), seci, thetabi, cfei);
				    }
				    phiCorFac1.add(Math.sqrt(min.fval()), seci, thetabi, 6);
				}
			}

			IndexedList<GraphErrors> histerrGroups_dtheta_vs_phi = new IndexedList<GraphErrors>(2);
			
			for(int ihisterrsectordth = 0; ihisterrsectordth < 6; ihisterrsectordth++) {
				for(int itheta_bindth = 0; itheta_bindth < 10; itheta_bindth++) {
					GraphErrors hedth = new GraphErrors();
					histerrGroups_dtheta_vs_phi.add(hedth, ihisterrsectordth, itheta_bindth);
				}
			}

			for(int seci = 0; seci < 6; seci++)
			{
				for(int thetadthi = 0; thetadthi < 9; thetadthi++)
				{
					for(int phidthi = 0; phidthi < 4; phidthi++)
					{
						histerrGroups_dtheta_vs_phi.getItem(seci, thetadthi).addPoint(phirotba[phidthi], mean1.getItem(seci, thetadthi,  phidthi),
																						0, sigma1.getItem(seci, thetadthi, phidthi));
					}
				}
			}
			
			for(int framehgedthi = 0; framehgedthi < 9; framehgedthi++){
				JFrame framedthphibin = new JFrame("#delta#theta vs. #phi Binned");
				framedthphibin.setSize(1500, 1000);
				EmbeddedCanvas canhgedth = new EmbeddedCanvas();
				framedthphibin.add(canhgedth);
				framedthphibin.setLocationRelativeTo(null);
				framedthphibin.setVisible(true);
				canhgedth.divide(3, 2);
				for(int seci = 0; seci < 6; seci++)
				{
					canhgedth.cd(seci);
					canhgedth.setFont("Arial");
					histerrGroups_dtheta_vs_phi.getItem(seci, framehgedthi).setTitle("Sector " + (seci+1) + " e^-: "
																						+ theta_bnd[framehgedthi]
																						+ " #degree < #theta <= "
																						+ theta_bnd[framehgedthi+1] + " #degree, "
																						+ " #delta#theta vs. #phi");
					histerrGroups_dtheta_vs_phi.getItem(seci, framehgedthi).setTitleX("#phi_midsector [#degree]");
					histerrGroups_dtheta_vs_phi.getItem(seci, framehgedthi).setTitleY("#delta#theta");
					histerrGroups_dtheta_vs_phi.getItem(seci, framehgedthi).setMarkerSize(5);
					histerrGroups_dtheta_vs_phi.getItem(seci, framehgedthi).setMarkerStyle(1);
					histerrGroups_dtheta_vs_phi.getItem(seci, framehgedthi).setLineThickness(1);
					canhgedth.getPad(seci).setTitleFontSize(20);
					canhgedth.getPad(seci).setAxisTitleFontSize(20);
					canhgedth.getPad(seci).setAxisLabelFontSize(20);
					canhgedth.getPad(seci).setAxisRange(-30, 30, 0.95, 1.05);
					canhgedth.draw(histerrGroups_dtheta_vs_phi.getItem(seci, framehgedthi), "same");
					F1D fdth = new F1D("fdp", "[a0]+([a1]*x)+([a2]*x*x)", -30, +30);
					fdth.setParameter(0, phiCorFac1.getItem(seci, framehgedthi, 0));
					fdth.setParameter(1, phiCorFac1.getItem(seci, framehgedthi, 1));
					fdth.setParameter(2, phiCorFac1.getItem(seci, framehgedthi, 2));
					fdth.setLineColor(2);
					fdth.setLineWidth(3);
					fdth.setOptStat(1110);
					canhgedth.draw(fdth, "same");
				}
			}
			
			IndexedList<Double> thetaCorFac1 = new IndexedList<>(3);
			
			for(int seci = 0; seci < 6; seci++)
			{
				for(int phcfi = 0; phcfi < 3; phcfi++)
				{
					List<Double> thetabcfi = new ArrayList<>();
					List<Double> cf1i = new ArrayList<>();
					List<Double> cfe1i = new ArrayList<>();
					
					for(int thetabfi = 0; thetabfi < 9; thetabfi++)
					{
						thetabcfi.add(thetamean[thetabfi]);
						cf1i.add(phiCorFac1.getItem(seci, thetabfi, phcfi));
						cfe1i.add(phiCorFac1.getItem(seci, thetabfi, phcfi+1));
					}
					
					FCN3 theFCN3 = new FCN3(thetabcfi, cf1i, cfe1i);
						
					MnUserParameters upar1 = new MnUserParameters();
					upar1.add("b0", 0, 0.01);
					upar1.add("b1", 0, 0.001);
					upar1.add("b2", 0, 0.0001);
					upar1.add("b3", 0, 0.00001);
					upar1.add("b4", 0, 0.000001);
					System.out.println("Initial parameters: "+upar1); 
				    
				    System.out.println("start migrad");
					MnMigrad migrad1 = new MnMigrad(theFCN3, upar1); 
				    FunctionMinimum min1 = migrad1.minimize(); 
				    if(!min1.isValid()) 
				    { 
				    	//try with higher strategy 
					    	System.out.println("FM is invalid, try with strategy = 2."); 
					    	MnMigrad migrad21 = new MnMigrad(theFCN3, min1.userState(), new MnStrategy(2)); 
					    	min1 = migrad21.minimize(); 
				    }
				    System.out.println("minimum: "+min1);
				    for(int thcfi = 0; thcfi < 5; thcfi++)
				    {
				    	thetaCorFac1.add(min1.userParameters().value(thcfi), seci, phcfi, thcfi);
				    }
				    for(int thcfei = 5; thcfei < 10; thcfei++)
				    {
				    	thetaCorFac1.add(min1.userParameters().error(thcfei-5), seci, phcfi, thcfei);
				    }
				    thetaCorFac1.add(Math.sqrt(min1.fval()), seci, phcfi, 10);
				}
			}

			IndexedList<GraphErrors> histerrGroupsCorFac_theta = new IndexedList<GraphErrors>(2);
			
			for(int ihisterrsectorcorfac1 = 0; ihisterrsectorcorfac1 < 6; ihisterrsectorcorfac1++) {
				for(int itheta_bincorfac1 = 0; itheta_bincorfac1 < 5; itheta_bincorfac1++) {
					GraphErrors heWcorfac1 = new GraphErrors();
					histerrGroupsCorFac_theta.add(heWcorfac1, ihisterrsectorcorfac1, itheta_bincorfac1);
				}
			}
			
			for(int seci = 0; seci < 6; seci ++)
			{
				for(int cc1i = 0; cc1i < 3; cc1i++)
				{
					GraphErrors hhecf1 = histerrGroupsCorFac_theta.getItem(seci, cc1i);
					for(int thetabgci = 0; thetabgci < 9; thetabgci++)
					{
						hhecf1.addPoint(thetamean[thetabgci], phiCorFac1.getItem(seci, thetabgci, cc1i),
										0, phiCorFac1.getItem(seci, thetabgci, cc1i+1));
					}
				}
			}
			
			for(int framehgecf1i = 0; framehgecf1i < 3; framehgecf1i++){
				JFrame frameWphibincorfac1 = new JFrame("#delta#theta: #phi Correction Coefficients vs. #theta Binned");
				frameWphibincorfac1.setSize(1500, 1000);
				EmbeddedCanvas canhgecf1 = new EmbeddedCanvas();
				frameWphibincorfac1.add(canhgecf1);
				frameWphibincorfac1.setLocationRelativeTo(null);
				frameWphibincorfac1.setVisible(true);
				canhgecf1.divide(3, 2);
				for(int seci = 0; seci < 6; seci++)
				{
					canhgecf1.cd(seci);
					canhgecf1.setFont("Arial");
					histerrGroupsCorFac_theta.getItem(seci, framehgecf1i).setTitle("Sector " + (seci+1) + " e^- #delta#theta: #phi Correction Coefficient a" + framehgecf1i + " vs. #theta");
					histerrGroupsCorFac_theta.getItem(seci, framehgecf1i).setTitleY("a" + framehgecf1i + " [arb. unit]");
					histerrGroupsCorFac_theta.getItem(seci, framehgecf1i).setMarkerSize(5);
					histerrGroupsCorFac_theta.getItem(seci, framehgecf1i).setMarkerStyle(1);
					histerrGroupsCorFac_theta.getItem(seci, framehgecf1i).setLineThickness(1);
					canhgecf1.getPad(seci).setTitleFontSize(20);
					canhgecf1.getPad(seci).setAxisTitleFontSize(20);
					canhgecf1.getPad(seci).setAxisLabelFontSize(20);
					if(framehgecf1i == 0) canhgecf1.getPad(seci).setAxisRange(7, 30, 0.9, 1.1);
					if(framehgecf1i == 1) canhgecf1.getPad(seci).setAxisRange(7, 30, -0.0025, 0.0025);
					if(framehgecf1i == 2) canhgecf1.getPad(seci).setAxisRange(7, 30, -0.005, 0.005);
					canhgecf1.draw(histerrGroupsCorFac_theta.getItem(seci, framehgecf1i), "same");
					F1D fcf1 = new F1D("fcf1", "[a0]+([a1]*x)+([a2]*x*x)+([a3]*x*x*x)+([a4]*x*x*x*x)", 7, 33);
					fcf1.setParameter(0, thetaCorFac1.getItem(seci, framehgecf1i, 0));
					fcf1.setParameter(1, thetaCorFac1.getItem(seci, framehgecf1i, 1));
					fcf1.setParameter(2, thetaCorFac1.getItem(seci, framehgecf1i, 2));
					fcf1.setParameter(3, thetaCorFac1.getItem(seci, framehgecf1i, 3));
					fcf1.setParameter(4, thetaCorFac1.getItem(seci, framehgecf1i, 4));
					fcf1.setLineColor(2);
					fcf1.setLineWidth(3);
					fcf1.setOptStat(111110);
					canhgecf1.draw(fcf1, "same");
				}
			}
			
			readerc1.open("C:/Users/joshtanj/Documents/download/merged_skim_e_bank_6535MeV_skim_elastic_5950_5899_5963_5958_5885.hipo");
			
			int corEventCounter1 = 0;
			while(readerc1.hasEvent())// && corEventCounter < 500000)
			{
				corEventCounter1++;
				processthetapCorrectedEvent(readerc1.getNextEvent(), corEventCounter, setsector, thetaCorFac, thetaCorFac1);
				if(corEventCounter1%50000 == 0) System.out.println("Corrected 1 Event: " + corEventCounter1);
			}
			
			readerc1.close();
			
			JFrame frameWcor1 = new JFrame("p, #theta Corrected W");
			frameWcor1.setSize(1500, 1000);
			EmbeddedCanvas canhg1cor1 = new EmbeddedCanvas();
			frameWcor1.add(canhg1cor1);
			frameWcor1.setLocationRelativeTo(null);
			frameWcor1.setVisible(true);
			canhg1cor1.divide(3, 2);
			for(int canhg1cor1i = 0; canhg1cor1i < 6; canhg1cor1i++){
				canhg1cor1.cd(canhg1cor1i);
				canhg1cor1.setFont("Arial");
				histGroups_theta_p_cor_W_sec.getItem(canhg1cor1i).setTitle("Sector " + (canhg1cor1i+1) + " W");
				histGroups_theta_p_cor_W_sec.getItem(canhg1cor1i).setTitleX("W [GeV]");
				histGroups_theta_p_cor_W_sec.getItem(canhg1cor1i).setTitleY("Counts");
				histGroups_theta_p_cor_W_sec.getItem(canhg1cor1i).setOptStat(10);
				canhg1cor1.getPad(canhg1cor1i).setTitleFontSize(20);
				canhg1cor1.getPad(canhg1cor1i).setAxisTitleFontSize(20);
				canhg1cor1.getPad(canhg1cor1i).setAxisLabelFontSize(20);
				canhg1cor1.getPad(canhg1cor1i).setStatBoxFontSize(20);
				canhg1cor1.draw(histGroups_theta_p_cor_W_sec.getItem(canhg1cor1i), "same");
				F1D f001c1 = new F1D("f001c1", "[amp]*gaus(x,[mean],[sigma])",
										0.7+(histGroups_theta_p_cor_W_sec.getItem(canhg1cor1i).getMaximumBin()*(1.2-0.7)/200)-0.030,
										0.7+(histGroups_theta_p_cor_W_sec.getItem(canhg1cor1i).getMaximumBin()*(1.2-0.7)/200)+0.020);
				f001c1.setParameter(0, histGroups_theta_p_cor_W_sec.getItem(canhg1cor1i).getMax()*1.2);
				f001c1.setParameter(1, 0.7+(histGroups_theta_p_cor_W_sec.getItem(canhg1cor1i).getMaximumBin()*(1.2-0.7)/200));
				f001c1.setParameter(2, 0.01);
				DataFitter.fit(f001c1, histGroups_theta_p_cor_W_sec.getItem(canhg1cor1i), "Q");
				f001c1.setLineColor(2);
				f001c1.setLineWidth(3);
				f001c1.setOptStat(11110);
				canhg1cor1.draw(f001c1, "same");
				Wmean.add(f001c1.getParameter(1), 2, canhg1cor1i);
				Wsigma.add(f001c1.getParameter(2), 2, canhg1cor1i);
			}
			
			JFrame framecor1Wsec1 = new JFrame("All Sector p, #theta Corrected W");
			framecor1Wsec1.setSize(500, 500);
			EmbeddedCanvas canhs1c1 = new EmbeddedCanvas();
			framecor1Wsec1.add(canhs1c1);
			framecor1Wsec1.setLocationRelativeTo(null);
			framecor1Wsec1.setVisible(true);
			canhs1c1.divide(1, 1);
			canhs1c1.cd(0);
			canhs1c1.setFont("Arial");
			hthpcorWs.setTitle("All Sector p, #theta Corrected W");
			hthpcorWs.setTitleX("W [GeV]");
			hthpcorWs.setTitleY("Counts");
			hthpcorWs.setOptStat(10);
			canhs1c1.getPad(0).setTitleFontSize(20);
			canhs1c1.getPad(0).setAxisTitleFontSize(20);
			canhs1c1.getPad(0).setAxisLabelFontSize(20);
			canhs1c1.getPad(0).setStatBoxFontSize(20);
			canhs1c1.draw(hthpcorWs, "same");
			F1D f00s1c1 = new F1D("f00s1c1", "[amp]*gaus(x,[mean],[sigma])",
									0.7+(hthpcorWs.getMaximumBin()*(1.2-0.7)/200)-0.030,
									0.7+(hthpcorWs.getMaximumBin()*(1.2-0.7)/200)+0.020);
			f00s1c1.setParameter(0, hthpcorWs.getMax()*1.2);
			f00s1c1.setParameter(1, 0.7+hthpcorWs.getMaximumBin()*(1.2-0.7)/200);
			f00s1c1.setParameter(2, 0.01);
			DataFitter.fit(f00s1c1, hthpcorWs, "Q");
			f00s1c1.setLineColor(2);
			f00s1c1.setLineWidth(3);
			f00s1c1.setOptStat(11110);
			canhs1c1.draw(f00s1c1, "same");
			
			IndexedList<Double> meancor1 = new IndexedList<>(3);
			IndexedList<Double> sigmacor1 = new IndexedList<>(3);
			
			for(int seci = 0; seci < 6; seci++)
			{
				for(int theta_binci = 0; theta_binci < 9; theta_binci++)
				{
					for(int phi_binci = 0; phi_binci < 4; phi_binci++)
					{
						F1D f003c1 = new F1D("f003c1", "[amp]*gaus(x,[mean],[sigma])",
												0.7+(histGroups_theta_p_cor_W_theta_phi_bin.getItem(seci, theta_binci, phi_binci).getMaximumBin()*(1.2-0.7)/200)-0.055,
												0.7+(histGroups_theta_p_cor_W_theta_phi_bin.getItem(seci, theta_binci, phi_binci).getMaximumBin()*(1.2-0.7)/200)+0.035);
						f003c1.setParameter(0, histGroups_theta_p_cor_W_theta_phi_bin.getItem(seci, theta_binci, phi_binci).getMax()*1.2);
						f003c1.setParameter(1, 0.7+histGroups_theta_p_cor_W_theta_phi_bin.getItem(seci, theta_binci, phi_binci).getMaximumBin()*(1.2-0.7)/200);
						f003c1.setParameter(2, 0.01);
						DataFitter.fit(f003c1, histGroups_theta_p_cor_W_theta_phi_bin.getItem(seci, theta_binci, phi_binci), "Q");
						System.out.println(f003c1.getParameter(1) + "   " + f003c1.getParameter(2));
						meancor1.add(f003c1.getParameter(1), seci, theta_binci, phi_binci);
						sigmacor1.add(f003c1.getParameter(2), seci, theta_binci, phi_binci);
					}
				}
			}
			
		System.out.println("p Correction Parameters");
		for(int seci = 0; seci < 6; seci++)
		{
			System.out.println("Sector " + (seci+1));
			for(int phcfi = 0; phcfi < 2; phcfi++)
			{
				System.out.println("B_" + phcfi + "_0, B_" + phcfi + "_1, B_" + phcfi + "_2, B_" + phcfi + "_3, B_" + phcfi + "_4");
			    System.out.println(thetaCorFac.getItem(seci, phcfi, 0) + ", "
				   					+ thetaCorFac.getItem(seci, phcfi, 1) + ", "
				   					+ thetaCorFac.getItem(seci, phcfi, 2) + ", "
				    				+ thetaCorFac.getItem(seci, phcfi, 3) + ", "
				    				+ thetaCorFac.getItem(seci, phcfi, 4));
			}
		}
		
		System.out.println("#theta Correction Parameters");
		for(int seci = 0; seci < 6; seci++)
		{
			for(int phcfi = 0; phcfi < 3; phcfi++)
			{
				System.out.println("b_" + phcfi + "_0, b_" + phcfi + "_1, b_" + phcfi + "_2, b_" + phcfi + "_3, b_" + phcfi + "_4");
				System.out.println(thetaCorFac1.getItem(setsector-1, phcfi, 0) + ", "
									+ thetaCorFac.getItem(seci, phcfi, 2) + ", "
									+ thetaCorFac.getItem(seci, phcfi, 3) + ", "
									+ thetaCorFac.getItem(seci, phcfi, 4) + ", "
									+ thetaCorFac.getItem(seci, phcfi, 5) + ", ");
			}
		}
	}
}