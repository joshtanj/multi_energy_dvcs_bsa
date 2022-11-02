package multi_energy_dvcs_bsa;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.swing.JFrame;

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
import org.slf4j.Marker;
import org.jlab.io.evio.EvioDataBank;
import org.jlab.io.hipo.HipoDataBank;

public class ecal_theta_phi_binner {
	
	static HipoDataSource reader = new HipoDataSource();
	static H2F hecal_det_y_vs_x = new H2F("ecal_det_y_vs_x", "ecal_det_y_vs_x", 500, -500, 500, 500, -500, 500);
	static H1F hecal_det_theta = new H1F("ecal_det_theta", "ecal_det_theta", 200, 3, 38);
	static H2F hecal_det_theta_vs_phi = new H2F("ecal_det_theta_vs_phi", "ecal_det_theta_vs_phi", 1500, -150, 210, 500, 3, 38);
	static H2F hecal_det_theta_vs_phi_rot = new H2F("ecal_det_theta_vs_phi_rot", "ecal_det_theta_vs_phi_rot", 500, -60, 60, 500, 3, 38);
	static H2F hecal_det_theta_vs_phi_rot_mir = new H2F("ecal_det_theta_vs_phi_rot_mir", "ecal_det_theta_vs_phi_rot_mir", 500, -60, 60, 500, 3, 38);
	static H1F hecal_det_phi_rot_mir_prj = new H1F("ecal_det_phi_rot_mir_prj", "ecal_det_phi_rot_mir_prj", 200, -30, 30);
	
	static IndexedList<H1F> histGroups_ecal_det_theta_bin= new IndexedList<H1F>(1);
	public static void ecal_det_theta_bin_histos() {
		for(int ibin = 0; ibin < 9; ibin++) {
			H1F ecal_det_theta_bin = new H1F("ecal_det_theta_bin", "ecal_det_theta_bin", 240, 3, 38);
			histGroups_ecal_det_theta_bin.add(ecal_det_theta_bin, ibin);
		}
	}
	
	static H1F hecal_det_phi_rot_mir_prj_bin1 = new H1F("ecal_det_phi_rot_mir_prj_bin1", "ecal_det_phi_rot_mir_prj_bin1", 200, -30, 30);
	static H1F hecal_det_phi_rot_mir_prj_bin2 = new H1F("ecal_det_phi_rot_mir_prj_bin2", "ecal_det_phi_rot_mir_prj_bin2", 200, -30, 30);
	static H1F hecal_det_phi_rot_mir_prj_bin3 = new H1F("ecal_det_phi_rot_mir_prj_bin3", "ecal_det_phi_rot_mir_prj_bin3", 200, -30, 30);
	
	static IndexedList<H1F> histGroups_ecal_det_phi_rot_mir_prj_theta_bin= new IndexedList<H1F>(1);
	static IndexedList<H1F> histGroups_ecal_det_phi_rot_mir_prj_bin= new IndexedList<H1F>(2);	
	
	public static void ecal_det_phi_rot_mir_prj_bin_histos() {
		for(int ithetabin = 0; ithetabin < 9; ithetabin++) {
			H1F ecal_det_phi_rot_mir_prj_theta_bin = new H1F("ecal_det_phi_rot_mir_prj_theta_bin", "ecal_det_phi_rot_mir_prj_theta_bin", 200, -30, 30);
			histGroups_ecal_det_phi_rot_mir_prj_theta_bin.add(ecal_det_phi_rot_mir_prj_theta_bin, ithetabin);
			for(int iphibin = 0; iphibin < 3; iphibin++) {
				H1F ecal_det_phi_rot_mir_prj_bin = new H1F("ecal_det_phi_rot_mir_prj_bin", "ecal_det_phi_rot_mir_prj_bin", 200, -30, 30);
				histGroups_ecal_det_phi_rot_mir_prj_bin.add(ecal_det_phi_rot_mir_prj_bin, ithetabin, iphibin);
			}
		}
	}
	
	static double dtr = Math.PI/180;
	
	static void processEvent(DataEvent event, int eventCounter)
	{
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Track") && event.hasBank("REC::Traj"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int ecount = 0;
			int protoncount = 0;
			int e_pindex = -1;
			int proton_pindex = -1;
			double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0.000511);
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0.938272);
			double E = 6.535;
		//	double E = 7.546;
			Vector3D v_e = new Vector3D (0, 0, 0);
			float theta_deg_e = 0;
			Vector3D v_proton = new Vector3D (0, 0, 0);
			float cdchi2_proton = 100;
			short cdNDF_proton = 100;
			int gammacount = 0;
			int gamma_pindex = -1;		
			LorentzVector p_gamma = new LorentzVector(0, 0, 0, 0);
			float theta_deg_gamma = 0;
			float phi_deg_gamma = -200;
			float phi_rot_deg_gamma = -200;
			double cX = 0;
			double cY = 0;
			boolean fid_ecal_gamma = false;
			double[] cX_split  = new double[] {87, 82, 85, 77, 78, 82};
			double[] t_left  = new double[] {58.7356, 62.8204, 62.2296, 53.7756, 58.2888, 54.5822};
			double[] t_right  = new double[] {58.7477, 51.2589, 59.2357, 56.2415, 60.8219, 49.8914};
			double[] s_left  = new double[] {0.582053, 0.544976, 0.549788, 0.56899, 0.56414, 0.57343};
			double[] s_right  = new double[] {-0.591876, -0.562926, -0.562246, -0.563726, -0.568902, -0.550729};
			double[] r_left  = new double[] {64.9348, 64.7541, 67.832, 55.9324, 55.9225, 60.0997};
			double[] r_right  = new double[] {65.424, 54.6992, 63.6628, 57.8931, 56.5367, 56.4641};
			double[] q_left  = new double[] {0.745578, 0.606081, 0.729202, 0.627239, 0.503674, 0.717899};
			double[] q_right  = new double[] {-0.775022, -0.633863, -0.678901, -0.612458, -0.455319, -0.692481};
			Vector3D dc_hit_e = new Vector3D (0, 0, 0);
			Vector3D ecal_hit_gamma = new Vector3D (0, 0, 0);
			Vector3D ecal_hit_rot_gamma = new Vector3D (0, 0, 0);
			for(int i = 0; i < rec.rows(); i++)
			{
				int pid = rec.getInt("pid", i);
				float px = rec.getFloat("px", i);
				float py = rec.getFloat("py", i);
				float pz = rec.getFloat("pz", i);
				float vx = rec.getFloat("vx", i);
				float vy = rec.getFloat("vy", i);
				float vz = rec.getFloat("vz", i);
				byte charge = rec.getByte("charge", i);
				float beta = rec.getFloat("beta", i);
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
				if(pid == 22)
				{
					gammacount++;
					p_gamma = new LorentzVector(px, py, pz, p);
					gamma_pindex = i;
					
					//Start: gamma fiducial cut
					if(event.hasBank("REC::Calorimeter"))
					{
						HipoDataBank reccal = (HipoDataBank) event.getBank("REC::Calorimeter");
						for(int l = 0; l < reccal.rows(); l++)
						{
							short pindex = reccal.getShort("pindex", l);
							byte detector = reccal.getByte("detector", l);
							byte layer = reccal.getByte("layer", l);
							float x = reccal.getFloat("x", l);
							float y = reccal.getFloat("y", l);
							float z = reccal.getFloat("z", l);
							float lu = reccal.getFloat("lu", l);
							float lv = reccal.getFloat("lv", l);
							float lw = reccal.getFloat("lw", l);
							byte sector = reccal.getByte("sector", l);
							if(pindex == gamma_pindex && detector == 7 && layer == 1)
							{
								ecal_hit_gamma = new Vector3D (x, y, z);
								theta_deg_gamma = (float) (57.3*ecal_hit_gamma.theta());
								phi_deg_gamma = (float) (57.3*ecal_hit_gamma.phi());
								if(phi_deg_gamma <= -150) phi_deg_gamma = phi_deg_gamma+360;
								ecal_hit_rot_gamma = new Vector3D (x*Math.cos(phi_rot[sector-1]/57.3)+y*Math.sin(phi_rot[sector-1]/57.3),
																y*Math.cos(phi_rot[sector-1]/57.3)-x*Math.sin(phi_rot[sector-1]/57.3), z);
								phi_rot_deg_gamma = (float) (57.3*ecal_hit_rot_gamma.phi());
							}
							if(pindex == gamma_pindex && detector == 7 && layer == 1)
							{
								cX = x*Math.cos(dtr*phi_rot[sector-1])+y*Math.sin(dtr*phi_rot[sector-1]);
								cY = y*Math.cos(dtr*phi_rot[sector-1])-x*Math.sin(dtr*phi_rot[sector-1]);
								if(cX > cX_split[sector-1] && cY < s_left[sector-1]*(cX-t_left[sector-1])
									&& cY > s_right[sector-1]*(cX-t_right[sector-1]))
									fid_ecal_gamma = true;
								else if(cX < cX_split[sector-1] && cY < q_left[sector-1]*(cX-r_left[sector-1])
											&& cY > q_right[sector-1]*(cX-r_right[sector-1]))
									fid_ecal_gamma = true;
							}
							
							if(pindex == gamma_pindex && detector == 7 && sector == 1)
							{
								if(layer == 1)
								{
									if((lw > 74.5 && lw < 79.5) || (lw > 85.5 && lw < 90.5) || (lw > 213 && lw < 218)
											|| (lw > 224.5 && lw < 229.5) || (lw > 410)) fid_ecal_gamma = false;
								}
								else if(layer == 4)
								{
									if((lv > 70 && lv < 93)) fid_ecal_gamma = false;
								}
								else if(layer == 7)
								{
									if((lu > 410.5)) fid_ecal_gamma = false;
									if((lv < 24) || (lv > 36.5 && lv < 60) || (lv > 411)) fid_ecal_gamma = false;
									if((lw < 21.5)) fid_ecal_gamma = false;
								}
							}
							else if(pindex == gamma_pindex && detector == 7 && sector == 2)
							{
								if(layer == 1)
								{
									if((lv > 102 && lv < 113)) fid_ecal_gamma = false;
								}
								else if(layer == 4)
								{
									if((lu > 396)) fid_ecal_gamma = false;
									if((lw > 363)) fid_ecal_gamma = false;
								}
								else if(layer == 7)
								{
									if((lu < 12)) fid_ecal_gamma = false;
									if((lw < 10.5) || (lw > 376)) fid_ecal_gamma = false;
								}
							}
							else if(pindex == gamma_pindex && detector == 7 && sector == 3)
							{
								if(layer == 1)
								{
									if((lv > 354.5 && lv < 376.5)) fid_ecal_gamma = false;
								}
								else if(layer == 4)
								{
									if((lu < 23)) fid_ecal_gamma = false;
									if((lw < 10) || (lw > 363)) fid_ecal_gamma = false;
								}
								else if(layer == 7)
								{
									if((lw > 387)) fid_ecal_gamma = false;
								}
							}
							else if(pindex == gamma_pindex && detector == 7 && sector == 4)
							{
								if(layer == 1)
								{
									if((lv < 13) || (lv > 230 && lv < 240.5)) fid_ecal_gamma = false;
									if((lw > 410)) fid_ecal_gamma = false;
								}
								else if(layer == 4)
								{
									if((lu < 20.5)) fid_ecal_gamma = false;
								}
								else if(layer == 7)
								{
									if((lw < 32.5)) fid_ecal_gamma = false;
								}
							}
							else if(pindex == gamma_pindex && detector == 7 && sector == 5)
							{
								if(layer == 4)
								{
									if((lv < 23)) fid_ecal_gamma = false;
									if((lw < 10)) fid_ecal_gamma = false;
								}
								else if(layer == 7)
								{
									if((lu > 193.5 && lu < 217)) fid_ecal_gamma = false;
									if((lv < 24)) fid_ecal_gamma = false;
								}
							}
							else if(pindex == gamma_pindex && detector == 7 && sector == 6)
							{
								if(layer == 1)
								{
									if((lw > 174.5 && lw < 179.5) || (lw > 185.5 && lw < 190.5)) fid_ecal_gamma = false;
								}
								else if(layer == 4)
								{
									if((lv < 11.5)) fid_ecal_gamma = false;
									if((lu < 20.5)) fid_ecal_gamma = false;
								}
								else if(layer == 7)
								{
									if((lv < 12) || (lv > 423)) fid_ecal_gamma = false;
									if((lw < 32.5)) fid_ecal_gamma = false;
								}
							}		
						}
					}
					// End: gamma fiducial cut
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
				}
			}
			
			LorentzVector p_B = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			LorentzVector p_T = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector q = new LorentzVector(0, 0, 0, 0);
			q.add(p_B);
			q.sub(p_e);
			double Q2 = -1*q.mass2();
			LorentzVector w = new LorentzVector(0, 0, 0, 0);
			w.add(p_B);
			w.add(p_T);
			w.sub(p_e);
			double W = w.mass();
			
			if(ecount == 1 && protoncount == 1 && gammacount == 1
					&& v_e.z() > -12 && v_e.z() < 7 && Math.abs(v_e.z()-v_proton.z()) < 2.5 + (2.5/p_proton.p())
					&& p_gamma.vect().theta(p_e.vect()) > 5
					&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
					&& 57.3*p_proton.theta() < 75
					&& fid_ecal_gamma == true
					&& theta_deg_e > 7 && theta_deg_e <= 30)
			{
				hecal_det_y_vs_x.fill(ecal_hit_gamma.x(), ecal_hit_gamma.y());
				hecal_det_theta.fill(theta_deg_gamma);
				hecal_det_theta_vs_phi.fill(phi_deg_gamma, theta_deg_gamma);
				hecal_det_theta_vs_phi_rot.fill(phi_rot_deg_gamma, theta_deg_gamma);
				hecal_det_theta_vs_phi_rot_mir.fill(Math.abs(phi_rot_deg_gamma), theta_deg_gamma);
				hecal_det_phi_rot_mir_prj.fill(Math.abs(phi_rot_deg_gamma));
				double[] theta_bnd  = new double[] {5, 7, 9, 11, 14, 17, 20, 25, 30, 35};
				for(int thbin = 0; thbin < 9; thbin++)
				{
					if(theta_deg_gamma > theta_bnd[thbin] && theta_deg_gamma <= theta_bnd[thbin+1])
					{
						histGroups_ecal_det_theta_bin.getItem(thbin).fill(theta_deg_gamma);
						histGroups_ecal_det_phi_rot_mir_prj_theta_bin.getItem(thbin).fill(Math.abs(phi_rot_deg_gamma));
						break;
					}
				}
			}
		}
	}
	
	public static void main(String[] args) {
		
		ecal_det_theta_bin_histos();
		ecal_det_phi_rot_mir_prj_bin_histos();
	
	//	reader.open("C:/Users/joshtanj/Documents/download/merged_skim_e_bank_6535MeV_skim_elastic_5950_5899_5963_5958_5885.hipo");
		reader.open("C:/Users/joshtanj/Documents/download/skim_6bank_epg_merged_6535MeV_skim16.hipo");
		
		int eventCounter = 1;
		while(reader.hasEvent())// && eventCounter < 1000000)
		{
			processEvent(reader.getNextEvent(), eventCounter);
			if(eventCounter%50000 == 0) System.out.println("Event: " + eventCounter);
			eventCounter++;
		}
		int Nevent = eventCounter-1;
		System.out.println("Number of events: " + Nevent);
		
		JFrame frame_ecal_hit = new JFrame("ECAL Hits");
		frame_ecal_hit.setSize(500, 500);
		EmbeddedCanvas can_ecal_hit = new EmbeddedCanvas();
		frame_ecal_hit.add(can_ecal_hit);
		frame_ecal_hit.setLocationRelativeTo(null);
		frame_ecal_hit.setVisible(true);
		can_ecal_hit.divide(1, 1);
		can_ecal_hit.cd(0);
		can_ecal_hit.setFont("Arial");
		hecal_det_y_vs_x.setTitle("ECAL y vs. x");
		hecal_det_y_vs_x.setTitleX("x [cm]");
		hecal_det_y_vs_x.setTitleY("y [cm]");
		can_ecal_hit.getPad(0).setTitleFontSize(32);
		can_ecal_hit.getPad(0).setAxisTitleFontSize(32);
		can_ecal_hit.getPad(0).setAxisLabelFontSize(24);
		can_ecal_hit.getPad(0).setStatBoxFontSize(18);
		can_ecal_hit.draw(hecal_det_y_vs_x, "same");
		
		JFrame frame_ecal_hit_theta = new JFrame("ECAL Hits #theta");
		frame_ecal_hit_theta.setSize(500, 500);
		EmbeddedCanvas can_ecal_hit_theta = new EmbeddedCanvas();
		frame_ecal_hit_theta.add(can_ecal_hit_theta);
		frame_ecal_hit_theta.setLocationRelativeTo(null);
		frame_ecal_hit_theta.setVisible(true);
		can_ecal_hit_theta.divide(1, 1);
		can_ecal_hit_theta.cd(0);
		can_ecal_hit_theta.setFont("Arial");
		hecal_det_theta.setTitle("ECAL #theta");
		hecal_det_theta.setTitleX("#theta [#degree]");
		hecal_det_theta.setTitleY("Counts");
		hecal_det_theta.setOptStat(10);
		can_ecal_hit_theta.getPad(0).setTitleFontSize(32);
		can_ecal_hit_theta.getPad(0).setAxisTitleFontSize(32);
		can_ecal_hit_theta.getPad(0).setAxisLabelFontSize(24);
		can_ecal_hit_theta.getPad(0).setStatBoxFontSize(18);
		can_ecal_hit_theta.draw(hecal_det_theta, "same");
		
		JFrame frame_ecal_hit_theta_vs_phi = new JFrame("ECAL Hits #theta vs. #phi");
		frame_ecal_hit_theta_vs_phi.setSize(1500, 500);
		EmbeddedCanvas can_ecal_hit_theta_vs_phi = new EmbeddedCanvas();
		frame_ecal_hit_theta_vs_phi.add(can_ecal_hit_theta_vs_phi);
		frame_ecal_hit_theta_vs_phi.setLocationRelativeTo(null);
		frame_ecal_hit_theta_vs_phi.setVisible(true);
		can_ecal_hit_theta_vs_phi.divide(1, 1);
		can_ecal_hit_theta_vs_phi.cd(0);
		can_ecal_hit_theta_vs_phi.setFont("Arial");
		hecal_det_theta_vs_phi.setTitle("ECAL #theta vs #phi");
		hecal_det_theta_vs_phi.setTitleX("#phi [#degree]");
		hecal_det_theta_vs_phi.setTitleY("#theta [#degree]");
		can_ecal_hit_theta_vs_phi.getPad(0).setTitleFontSize(32);
		can_ecal_hit_theta_vs_phi.getPad(0).setAxisTitleFontSize(32);
		can_ecal_hit_theta_vs_phi.getPad(0).setAxisLabelFontSize(24);
		can_ecal_hit_theta_vs_phi.getPad(0).setStatBoxFontSize(18);
		can_ecal_hit_theta_vs_phi.draw(hecal_det_theta_vs_phi, "same");
		
		JFrame frame_ecal_hit_theta_vs_phi_prj = new JFrame("ECAL Hits #phi Projections");
		frame_ecal_hit_theta_vs_phi_prj.setSize(1500, 500);
		EmbeddedCanvas can_ecal_hit_theta_vs_phi_prj = new EmbeddedCanvas();
		frame_ecal_hit_theta_vs_phi_prj.add(can_ecal_hit_theta_vs_phi_prj);
		frame_ecal_hit_theta_vs_phi_prj.setLocationRelativeTo(null);
		frame_ecal_hit_theta_vs_phi_prj.setVisible(true);
		can_ecal_hit_theta_vs_phi_prj.divide(3, 1);
		can_ecal_hit_theta_vs_phi_prj.cd(0);
		can_ecal_hit_theta_vs_phi_prj.setFont("Arial");
		hecal_det_theta_vs_phi_rot.setTitle("Superpositioned ECAL #theta vs. #phi");
		hecal_det_theta_vs_phi_rot.setTitleX("#phi [#degree]");
		hecal_det_theta_vs_phi_rot.setTitleY("#theta [#degree]");
		can_ecal_hit_theta_vs_phi_prj.getPad(0).setTitleFontSize(32);
		can_ecal_hit_theta_vs_phi_prj.getPad(0).setAxisTitleFontSize(32);
		can_ecal_hit_theta_vs_phi_prj.getPad(0).setAxisLabelFontSize(24);
		can_ecal_hit_theta_vs_phi_prj.getPad(0).setStatBoxFontSize(18);
		can_ecal_hit_theta_vs_phi_prj.draw(hecal_det_theta_vs_phi_rot, "same");
		can_ecal_hit_theta_vs_phi_prj.cd(1);
		can_ecal_hit_theta_vs_phi_prj.setFont("Arial");
		hecal_det_theta_vs_phi_rot_mir.setTitle("Superpositioned Mirrored ECAL #theta vs. #phi");
		hecal_det_theta_vs_phi_rot_mir.setTitleX("#phi [#degree]");
		hecal_det_theta_vs_phi_rot_mir.setTitleY("#theta [#degree]");
		can_ecal_hit_theta_vs_phi_prj.getPad(1).setTitleFontSize(32);
		can_ecal_hit_theta_vs_phi_prj.getPad(1).setAxisTitleFontSize(32);
		can_ecal_hit_theta_vs_phi_prj.getPad(1).setAxisLabelFontSize(24);
		can_ecal_hit_theta_vs_phi_prj.getPad(1).setStatBoxFontSize(18);
		can_ecal_hit_theta_vs_phi_prj.draw(hecal_det_theta_vs_phi_rot_mir, "same");
		can_ecal_hit_theta_vs_phi_prj.cd(2);
		can_ecal_hit_theta_vs_phi_prj.setFont("Arial");
		hecal_det_phi_rot_mir_prj.setTitle("Superpositioned Mirrored ECAL #phi");
		hecal_det_phi_rot_mir_prj.setTitleX("#phi [#degree]");
		hecal_det_phi_rot_mir_prj.setTitleY("Counts");
		hecal_det_phi_rot_mir_prj.setOptStat(10);
		can_ecal_hit_theta_vs_phi_prj.getPad(2).setTitleFontSize(32);
		can_ecal_hit_theta_vs_phi_prj.getPad(2).setAxisTitleFontSize(32);
		can_ecal_hit_theta_vs_phi_prj.getPad(2).setAxisLabelFontSize(24);
		can_ecal_hit_theta_vs_phi_prj.getPad(2).setStatBoxFontSize(18);
		can_ecal_hit_theta_vs_phi_prj.draw(hecal_det_phi_rot_mir_prj, "same");
		
		double[] theta_bnd  = new double[] {5, 7, 9, 11, 14, 17, 20, 25, 30, 35};
		
		JFrame frame_ecal_hit_theta_bin = new JFrame("ECAL Hits #theta Bins");
		frame_ecal_hit_theta_bin.setSize(1500, 1000);
		EmbeddedCanvas can_ecal_hit_theta_bin = new EmbeddedCanvas();
		frame_ecal_hit_theta_bin.add(can_ecal_hit_theta_bin);
		frame_ecal_hit_theta_bin.setLocationRelativeTo(null);
		frame_ecal_hit_theta_bin.setVisible(true);
		can_ecal_hit_theta_bin.divide(3, 3);	
		for(int thbin = 0; thbin < 9; thbin++)
		{
			can_ecal_hit_theta_bin.cd(thbin);
			can_ecal_hit_theta_bin.setFont("Arial");
			histGroups_ecal_det_theta_bin.getItem(thbin).setTitle("ECAL " + theta_bnd[thbin] + "#degree < #theta <= "
																	+ theta_bnd[thbin+1] + "#degree");
			histGroups_ecal_det_theta_bin.getItem(thbin).setTitleX("#theta [#degree]");
			histGroups_ecal_det_theta_bin.getItem(thbin).setTitleY("Counts");
			histGroups_ecal_det_theta_bin.getItem(thbin).setOptStat(110);
			can_ecal_hit_theta_bin.getPad(thbin).setTitleFontSize(32);
			can_ecal_hit_theta_bin.getPad(thbin).setAxisTitleFontSize(32);
			can_ecal_hit_theta_bin.getPad(thbin).setAxisLabelFontSize(24);
			can_ecal_hit_theta_bin.getPad(thbin).setStatBoxFontSize(18);
			can_ecal_hit_theta_bin.draw(histGroups_ecal_det_theta_bin.getItem(thbin), "same");
		}
		
		int phidiv = 2;
		
		for(int bini = 0; bini < 200; bini++){
			if(hecal_det_phi_rot_mir_prj_bin1.getEntries() < hecal_det_phi_rot_mir_prj.getEntries()/phidiv){
				for(int filli = 0; filli < hecal_det_phi_rot_mir_prj.getBinContent(bini); filli++){
					hecal_det_phi_rot_mir_prj_bin1.fill(-30+0.3*bini);
				}
			}
			else if(hecal_det_phi_rot_mir_prj_bin2.getEntries() < hecal_det_phi_rot_mir_prj.getEntries()/phidiv){
				for(int filli = 0; filli < hecal_det_phi_rot_mir_prj.getBinContent(bini); filli++){
					hecal_det_phi_rot_mir_prj_bin2.fill(-30+0.3*bini);
				}
			}
			else{
				for(int filli = 0; filli < hecal_det_phi_rot_mir_prj.getBinContent(bini); filli++){
					hecal_det_phi_rot_mir_prj_bin3.fill(-30+0.3*bini);
				}
			}
		}
		
		JFrame frame_ecal_hit_phi_bin = new JFrame("ECAL Hits #phi Bins");
		frame_ecal_hit_phi_bin.setSize(1500, 500);
		EmbeddedCanvas can_ecal_hit_phi_bin = new EmbeddedCanvas();
		frame_ecal_hit_phi_bin.add(can_ecal_hit_phi_bin);
		frame_ecal_hit_phi_bin.setLocationRelativeTo(null);
		frame_ecal_hit_phi_bin.setVisible(true);
		can_ecal_hit_phi_bin.divide(3, 1);
		can_ecal_hit_phi_bin.cd(0);
		can_ecal_hit_phi_bin.setFont("Arial");
		hecal_det_phi_rot_mir_prj_bin1.setTitle("ECAL #phi Bin 1");
		hecal_det_phi_rot_mir_prj_bin1.setTitleX("#phi [#degree]");
		hecal_det_phi_rot_mir_prj_bin1.setTitleY("Counts");
		hecal_det_phi_rot_mir_prj_bin1.setOptStat(110);
		can_ecal_hit_phi_bin.getPad(0).setTitleFontSize(32);
		can_ecal_hit_phi_bin.getPad(0).setAxisTitleFontSize(32);
		can_ecal_hit_phi_bin.getPad(0).setAxisLabelFontSize(24);
		can_ecal_hit_phi_bin.getPad(0).setStatBoxFontSize(18);
		can_ecal_hit_phi_bin.draw(hecal_det_phi_rot_mir_prj_bin1, "same");
		can_ecal_hit_phi_bin.cd(1);
		can_ecal_hit_phi_bin.setFont("Arial");
		hecal_det_phi_rot_mir_prj_bin2.setTitle("ECAL #phi Bin 2");
		hecal_det_phi_rot_mir_prj_bin2.setTitleX("#phi [#degree]");
		hecal_det_phi_rot_mir_prj_bin2.setTitleY("Counts");
		hecal_det_phi_rot_mir_prj_bin2.setOptStat(110);
		can_ecal_hit_phi_bin.getPad(1).setTitleFontSize(32);
		can_ecal_hit_phi_bin.getPad(1).setAxisTitleFontSize(32);
		can_ecal_hit_phi_bin.getPad(1).setAxisLabelFontSize(24);
		can_ecal_hit_phi_bin.getPad(1).setStatBoxFontSize(18);
		can_ecal_hit_phi_bin.draw(hecal_det_phi_rot_mir_prj_bin2, "same");
		can_ecal_hit_phi_bin.cd(2);
		can_ecal_hit_phi_bin.setFont("Arial");
		hecal_det_phi_rot_mir_prj_bin3.setTitle("ECAL #phi Bin 3");
		hecal_det_phi_rot_mir_prj_bin3.setTitleX("#phi [#degree]");
		hecal_det_phi_rot_mir_prj_bin3.setTitleY("Counts");
		hecal_det_phi_rot_mir_prj_bin3.setOptStat(110);
		can_ecal_hit_phi_bin.getPad(2).setTitleFontSize(32);
		can_ecal_hit_phi_bin.getPad(2).setAxisTitleFontSize(32);
		can_ecal_hit_phi_bin.getPad(2).setAxisLabelFontSize(24);
		can_ecal_hit_phi_bin.getPad(2).setStatBoxFontSize(18);
		can_ecal_hit_phi_bin.draw(hecal_det_phi_rot_mir_prj_bin3, "same");
		
		for(int thetabini = 0; thetabini < 9; thetabini++){
			for(int bini = 0; bini < 200; bini++){
				if(histGroups_ecal_det_phi_rot_mir_prj_bin.getItem(thetabini, 0).getEntries()
						< histGroups_ecal_det_phi_rot_mir_prj_theta_bin.getItem(thetabini).getEntries()/phidiv){
					for(int filli = 0; filli < histGroups_ecal_det_phi_rot_mir_prj_theta_bin.getItem(thetabini).getBinContent(bini); filli++){
						histGroups_ecal_det_phi_rot_mir_prj_bin.getItem(thetabini, 0).fill(-30+0.3*bini);
					}
				}
				else if(histGroups_ecal_det_phi_rot_mir_prj_bin.getItem(thetabini, 1).getEntries()
						< histGroups_ecal_det_phi_rot_mir_prj_theta_bin.getItem(thetabini).getEntries()/phidiv){
					for(int filli = 0; filli < histGroups_ecal_det_phi_rot_mir_prj_theta_bin.getItem(thetabini).getBinContent(bini); filli++){
						histGroups_ecal_det_phi_rot_mir_prj_bin.getItem(thetabini, 1).fill(-30+0.3*bini);
					}
				}
				else{
					for(int filli = 0; filli < histGroups_ecal_det_phi_rot_mir_prj_theta_bin.getItem(thetabini).getBinContent(bini); filli++){
						histGroups_ecal_det_phi_rot_mir_prj_bin.getItem(thetabini, 2).fill(-30+0.3*bini);
					}
				}
			}
		}
		
		for(int thetabini = 0; thetabini < 9; thetabini++){
			JFrame frame_ecal_theta_phi_bin = new JFrame("ECAL Hits #theta Bin " + (thetabini+1));
			frame_ecal_theta_phi_bin.setSize(1500, 500);
			EmbeddedCanvas can_ecal_theta_phi_bin = new EmbeddedCanvas();
			frame_ecal_theta_phi_bin.add(can_ecal_theta_phi_bin);
			frame_ecal_theta_phi_bin.setLocationRelativeTo(null);
			frame_ecal_theta_phi_bin.setVisible(true);
			can_ecal_theta_phi_bin.divide(3, 1);
			for(int phibini = 0; phibini < 3; phibini++){
				can_ecal_theta_phi_bin.cd(phibini);
				can_ecal_theta_phi_bin.setFont("Arial");
				histGroups_ecal_det_phi_rot_mir_prj_bin.getItem(thetabini, phibini).setTitle("ECAL #theta Bin " + (thetabini+1)
																									+ " #phi Bin " + (phibini+1));
				histGroups_ecal_det_phi_rot_mir_prj_bin.getItem(thetabini, phibini).setTitleX("#phi [#degree]");
				histGroups_ecal_det_phi_rot_mir_prj_bin.getItem(thetabini, phibini).setTitleY("Counts");
				histGroups_ecal_det_phi_rot_mir_prj_bin.getItem(thetabini, phibini).setOptStat(110);
				can_ecal_theta_phi_bin.getPad(phibini).setTitleFontSize(32);
				can_ecal_theta_phi_bin.getPad(phibini).setAxisTitleFontSize(32);
				can_ecal_theta_phi_bin.getPad(phibini).setAxisLabelFontSize(24);
				can_ecal_theta_phi_bin.getPad(phibini).setStatBoxFontSize(18);
				can_ecal_theta_phi_bin.draw(histGroups_ecal_det_phi_rot_mir_prj_bin.getItem(thetabini, phibini), "same");
			}
		}
	}
}