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

public class ftcal_kinematic_mapping {
	
	static H2F hftcal_det_y_vs_x = new H2F("ftcal_det_y_vs_x", "ftcal_det_y_vs_x", 500, -17.5, 17.5, 500, -17.5, 17.5);
	static H2F hftcal_det_theta_vs_phi = new H2F("ftcal_det_theta_vs_phi", "ftcal_det_theta_vs_phi", 500, -180, 180, 500, 0, 7);
	static H2F hftcal_det_cryst_y_vs_x = new H2F("ftcal_det_cryst_y_vs_x", "ftcal_det_cryst_y_vs_x", 24, -18.36, 18.36, 24, -18.36, 18.36);
	static H2F hftcal_det_cryst_cy_vs_cx = new H2F("ftcal_det_cryst_cy_vs_cx", "ftcal_det_cryst_cy_vs_cx", 24, -0.5, 23.5, 24, -0.5, 23.5);
	
	static IndexedList<H1F> histGroups_theta_cone_gamma = new IndexedList<H1F>(2);
	public static void theta_cone_gamma_histos() {
		for(int xbin = 0; xbin  < 24; xbin++) {
			for(int ybin = 0; ybin < 24; ybin++)
			{
				H1F htheta_cone_gamma = new H1F("theta_cone_gamma", "theta_cone_gamma", 100, 0, 5);
				histGroups_theta_cone_gamma.add(htheta_cone_gamma, xbin, ybin);
			}
		}
	}
	
	static IndexedList<H1F> histGroups_X_eg_m = new IndexedList<H1F>(2);
	public static void X_eg_m_histos() {
		for(int xbin = 0; xbin  < 24; xbin++) {
			for(int ybin = 0; ybin < 24; ybin++)
			{
				H1F hX_eg_m = new H1F("X_eg_m", "X_eg_m", 100, 0, 3);
				histGroups_X_eg_m.add(hX_eg_m, xbin, ybin);
			}
		}
	}
	
	static IndexedList<H1F> histGroups_X_epg_pt = new IndexedList<H1F>(2);
	public static void X_epg_pt_histos() {
		for(int xbin = 0; xbin  < 24; xbin++) {
			for(int ybin = 0; ybin < 24; ybin++)
			{
				H1F hX_epg_pt = new H1F("X_epg_pt", "X_epg_pt", 100, 0, 0.8);
				histGroups_X_epg_pt.add(hX_epg_pt, xbin, ybin);
			}
		}
	}
	
	static IndexedList<H1F> histGroups_X_epg_E = new IndexedList<H1F>(2);
	public static void X_epg_E_histos() {
		for(int xbin = 0; xbin  < 24; xbin++) {
			for(int ybin = 0; ybin < 24; ybin++)
			{
				H1F hX_epg_E = new H1F("X_epg_E", "X_epg_E", 100, -2, 2);
				histGroups_X_epg_E.add(hX_epg_E, xbin, ybin);
			}
		}
	}
	
	static HipoDataSource reader = new HipoDataSource();
	
	static double rc = Math.PI/180;
	
	static void processEvent(DataEvent event, int eventCounter)
	{
		if(event.hasBank("REC::Particle"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int ecount = 0;
			int e_pindex = -1;
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0);
			Vector3D v_e = new Vector3D (0, 0, 0);
			int protoncount = 0;
			int proton_pindex = -1;
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0);
			Vector3D v_proton = new Vector3D (0, 0, 0);
			float cdchi2_proton = 100;
			short cdNDF_proton = 100;
			int gammacount = 0;
			int gamma_pindex = -1;		
			LorentzVector p_gamma = new LorentzVector(0, 0, 0, 0);
			boolean fid_ftcal_gamma = false;
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
				if(pid == 22)
				{
					gammacount++;
					p_gamma = new LorentzVector(px, py, pz, p);
					gamma_pindex = i;
				}	
			}
			HipoDataBank rectrac = (HipoDataBank) event.getBank("REC::Track");
			byte proton_detector = 0;
			for(int j = 0; j < rectrac.rows(); j++)
			{
				short pindex = rectrac.getShort("pindex", j);
				byte detector = rectrac.getByte("detector", j);
				float chi2 = rectrac.getFloat("chi2", j);
				short NDF = rectrac.getShort("NDF", j);
				if(pindex == proton_pindex)
				{
					proton_detector = detector;
					cdchi2_proton = chi2;
					cdNDF_proton = NDF;
				}
			}
			
			double E = 7.546;
			LorentzVector p_B = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			LorentzVector p_T = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector p_X_ep = new LorentzVector(0, 0, 0, 0);
			p_X_ep.add(p_B);
			p_X_ep.add(p_T);
			p_X_ep.sub(p_e);
			p_X_ep.sub(p_proton);
			LorentzVector p_X_eg = new LorentzVector(0, 0, 0, 0);
			p_X_eg.add(p_B);
			p_X_eg.add(p_T);
			p_X_eg.sub(p_e);
			p_X_eg.sub(p_gamma);
			LorentzVector p_X_epg = new LorentzVector(0, 0, 0, 0);
			p_X_epg.add(p_B);
			p_X_epg.add(p_T);
			p_X_epg.sub(p_e);
			p_X_epg.sub(p_proton);
			p_X_epg.sub(p_gamma);
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
					&& Q2 > 1 && W > 2)
			{
				if(event.hasBank("REC::ForwardTagger"))
				{
					HipoDataBank recft = (HipoDataBank) event.getBank("REC::ForwardTagger");
					for(int m = 0; m < recft.rows(); m++)
					{
						short pindex = recft.getShort("pindex", m);
						byte detector = recft.getByte("detector", m);
						byte layer = recft.getByte("layer", m);
						float x = recft.getFloat("x", m);
						float y = recft.getFloat("y", m);
						float z = recft.getFloat("z", m);
						Vector3D ftcal_det = new Vector3D (x, y, z);
						double xtrans = x + 18.36;
						double ytrans = y + 18.36;
						int cx = (int) (xtrans/1.53);
						int cy = (int) (ytrans/1.53);
						if(pindex == gamma_pindex && detector == 10 && layer == 1)
						{
							if(ftcal_det.rho() > 8.8 && !(ftcal_det.rho() > 15.5)
								&& (x+8.5)*(x+8.5)+(y-10)*(y-10) > 1.5*1.5 
								&& (x+6)*(x+6)/(2.3*2.3)+(y+13)*(y+13)/(1.6*1.6) > 1
								&& (x-3.8)*(x-3.8)+(y+6.8)*(y+6.8) > 1.7*1.7) fid_ftcal_gamma = true;
							if(cx == 7 && cy == 15) fid_ftcal_gamma = false;
						}
						if(pindex == gamma_pindex && detector == 10 && layer == 1 && fid_ftcal_gamma == true)
						{
							double ftcal_det_phi = 180*ftcal_det.phi()/Math.PI;
							double ftcal_det_theta = 180*ftcal_det.theta()/Math.PI;
							hftcal_det_y_vs_x.fill(x, y);
							hftcal_det_theta_vs_phi.fill(ftcal_det_phi, ftcal_det_theta);
							hftcal_det_cryst_y_vs_x.fill(x, y);
							hftcal_det_cryst_cy_vs_cx.fill((cx-0.5), (cy-0.5));
							histGroups_theta_cone_gamma.getItem(cx, cy).fill(p_gamma.vect().theta(p_X_ep.vect()));
							histGroups_X_eg_m.getItem(cx, cy).fill(p_X_eg.mass());
							histGroups_X_epg_pt.getItem(cx, cy).fill(Math.sqrt(p_X_epg.px()*p_X_epg.px()+p_X_epg.py()*p_X_epg.py()));
							histGroups_X_epg_E.getItem(cx, cy).fill(p_X_epg.e());
						}
					}
				}
			}
		}
	}
	
	public static void main(String[] args) {
		
		theta_cone_gamma_histos();
		X_eg_m_histos();
		X_epg_pt_histos();
		X_epg_E_histos();
		
		reader.open("C:/Users/joshtanj/Documents/download/skim_epg_bank_merged_7546MeV_skim16.hipo");
		
		int eventCounter = 1;
		while(reader.hasEvent())// && eventCounter < 20000000)
		{
			processEvent(reader.getNextEvent(), eventCounter);
			if(eventCounter%50000 == 0) System.out.println("Event: " + eventCounter);
			eventCounter++;
		}
		
		JFrame frame_ftcal_hit = new JFrame("FTCAL Hits");
		frame_ftcal_hit.setSize(1500, 750);
		EmbeddedCanvas can_ftcal_hit = new EmbeddedCanvas();
		frame_ftcal_hit.add(can_ftcal_hit);
		frame_ftcal_hit.setLocationRelativeTo(null);
		frame_ftcal_hit.setVisible(true);
		can_ftcal_hit.divide(2, 1);
		can_ftcal_hit.cd(0);
		can_ftcal_hit.setFont("Arial");
		hftcal_det_y_vs_x.setTitle("FTCAL y vs. x");
		hftcal_det_y_vs_x.setTitleX("x [cm]");
		hftcal_det_y_vs_x.setTitleY("y [cm]");
		can_ftcal_hit.getPad(0).setTitleFontSize(32);
		can_ftcal_hit.getPad(0).setAxisTitleFontSize(32);
		can_ftcal_hit.getPad(0).setAxisLabelFontSize(24);
		can_ftcal_hit.draw(hftcal_det_y_vs_x, "same");
		can_ftcal_hit.cd(1);
		can_ftcal_hit.setFont("Arial");
		hftcal_det_theta_vs_phi.setTitle("FTCAL #theta vs. #phi");
		hftcal_det_theta_vs_phi.setTitleX("#phi [deg]");
		hftcal_det_theta_vs_phi.setTitleY("#theta [deg]");
		can_ftcal_hit.getPad(1).setTitleFontSize(32);
		can_ftcal_hit.getPad(1).setAxisTitleFontSize(32);
		can_ftcal_hit.getPad(1).setAxisLabelFontSize(24);
		can_ftcal_hit.draw(hftcal_det_theta_vs_phi, "same");
		
		JFrame frame_ftcal_cryst_hit = new JFrame("FTCAL Crystal Hits");
		frame_ftcal_cryst_hit.setSize(1500, 750);
		EmbeddedCanvas can_ftcal_cryst_hit = new EmbeddedCanvas();
		frame_ftcal_cryst_hit.add(can_ftcal_cryst_hit);
		frame_ftcal_cryst_hit.setLocationRelativeTo(null);
		frame_ftcal_cryst_hit.setVisible(true);
		can_ftcal_cryst_hit.divide(2, 1);
		can_ftcal_cryst_hit.cd(0);
		can_ftcal_cryst_hit.setFont("Arial");
		hftcal_det_cryst_y_vs_x.setTitle("FTCAL Crystal y vs. x");
		hftcal_det_cryst_y_vs_x.setTitleX("x [cm]");
		hftcal_det_cryst_y_vs_x.setTitleY("y [cm]");
		can_ftcal_cryst_hit.getPad(0).setTitleFontSize(32);
		can_ftcal_cryst_hit.getPad(0).setAxisTitleFontSize(32);
		can_ftcal_cryst_hit.getPad(0).setAxisLabelFontSize(24);
		can_ftcal_cryst_hit.draw(hftcal_det_cryst_y_vs_x, "same");
		can_ftcal_cryst_hit.cd(1);
		can_ftcal_cryst_hit.setFont("Arial");
		hftcal_det_cryst_cy_vs_cx.setTitle("FTCAL Crystalized y vs. x");
		hftcal_det_cryst_cy_vs_cx.setTitleX("x [crystal unit]");
		hftcal_det_cryst_cy_vs_cx.setTitleY("y [crystal unit]");
		can_ftcal_cryst_hit.getPad(1).setTitleFontSize(32);
		can_ftcal_cryst_hit.getPad(1).setAxisTitleFontSize(32);
		can_ftcal_cryst_hit.getPad(1).setAxisLabelFontSize(24);
		can_ftcal_cryst_hit.draw(hftcal_det_cryst_cy_vs_cx, "same");
		
		IndexedList<Double> mean_tc = new IndexedList<>(2);
		IndexedList<Double> sigma_tc = new IndexedList<>(2);
		IndexedList<Double> mean_pm = new IndexedList<>(2);
		IndexedList<Double> sigma_pm = new IndexedList<>(2);
		IndexedList<Double> mean_Xpt = new IndexedList<>(2);
		IndexedList<Double> sigma_Xpt = new IndexedList<>(2);
		IndexedList<Double> mean_Xe = new IndexedList<>(2);
		IndexedList<Double> sigma_Xe = new IndexedList<>(2);
		
		for(int cx = 0; cx  < 24; cx++)
		{
			for(int cy = 0; cy < 24; cy++)
			{
				if(histGroups_theta_cone_gamma.getItem(cx, cy).getEntries() > 10)
				{
					F1D ftc = new F1D("ftc", "[amp]*gaus(x,[mean],[sigma])",
										histGroups_theta_cone_gamma.getItem(cx, cy).getMaximumBin()*0.05, 
										histGroups_theta_cone_gamma.getItem(cx, cy).getRMS()*3);
					ftc.setParameter(0, histGroups_theta_cone_gamma.getItem(cx, cy).getMax());
					ftc.setParLimits(0, histGroups_theta_cone_gamma.getItem(cx, cy).getMax()+1,
										histGroups_theta_cone_gamma.getItem(cx, cy).getMax()-1);
					ftc.setParameter(1, (histGroups_theta_cone_gamma.getItem(cx, cy).getMaximumBin()*0.05));
					ftc.setParLimits(1, ((histGroups_theta_cone_gamma.getItem(cx, cy).getMaximumBin()-1)*0.05),
										((histGroups_theta_cone_gamma.getItem(cx, cy).getMaximumBin()+1)*0.05));
					ftc.setParameter(2, histGroups_theta_cone_gamma.getItem(cx, cy).getRMS()/2);
					ftc.setParLimits(2, histGroups_theta_cone_gamma.getItem(cx, cy).getRMS()/4,
										histGroups_theta_cone_gamma.getItem(cx, cy).getRMS());
					DataFitter.fit(ftc, histGroups_theta_cone_gamma.getItem(cx, cy), "Q");
					mean_tc.add(ftc.getParameter(1), cx, cy);
					sigma_tc.add(Math.abs(ftc.getParameter(2)), cx, cy);
					
					F1D fpm = new F1D("fpm", "[amp]*gaus(x,[mean],[sigma])", (histGroups_X_eg_m.getItem(cx, cy).getMaximumBin()*0.03-1),
										(histGroups_X_eg_m.getItem(cx, cy).getMaximumBin()*0.03+1));
					fpm.setParameter(0, histGroups_X_eg_m.getItem(cx, cy).getMax());
					fpm.setParameter(1, histGroups_X_eg_m.getItem(cx, cy).getMaximumBin()*0.03);
					fpm.setParameter(2, histGroups_X_eg_m.getItem(cx, cy).getRMS()/2);
					fpm.setParLimits(2, histGroups_X_eg_m.getItem(cx, cy).getRMS()/4,
										histGroups_X_eg_m.getItem(cx, cy).getRMS());
					DataFitter.fit(fpm, histGroups_X_eg_m.getItem(cx, cy), "Q");
					mean_pm.add(fpm.getParameter(1), cx, cy);
					sigma_pm.add(Math.abs(fpm.getParameter(2)), cx, cy);
					
					F1D fXpt = new F1D("fXpt", "[amp]*gaus(x,[mean],[sigma])", histGroups_X_epg_pt.getItem(cx, cy).getMaximumBin()*0.008,
										histGroups_X_epg_pt.getItem(cx, cy).getRMS()*3);
					fXpt.setParameter(0, histGroups_X_epg_pt.getItem(cx, cy).getMax());
					fXpt.setParLimits(0, histGroups_X_epg_pt.getItem(cx, cy).getMax()+1,
										histGroups_X_epg_pt.getItem(cx, cy).getMax()-1);
					fXpt.setParameter(1, (histGroups_X_epg_pt.getItem(cx, cy).getMaximumBin()*0.008));
					fXpt.setParLimits(1, ((histGroups_X_epg_pt.getItem(cx, cy).getMaximumBin()-1)*0.008),
									((histGroups_X_epg_pt.getItem(cx, cy).getMaximumBin()+1)*0.008));
					fXpt.setParameter(2, histGroups_X_epg_pt.getItem(cx, cy).getRMS()/2);
					fXpt.setParLimits(2, histGroups_X_epg_pt.getItem(cx, cy).getRMS()/4,
										histGroups_X_epg_pt.getItem(cx, cy).getRMS());
					DataFitter.fit(fXpt, histGroups_X_epg_pt.getItem(cx, cy), "Q");
					mean_Xpt.add(fXpt.getParameter(1), cx, cy);
					sigma_Xpt.add(Math.abs(fXpt.getParameter(2)), cx, cy);
					
					F1D fXe = new F1D("fXe", "[amp]*gaus(x,[mean],[sigma])", (histGroups_X_epg_E.getItem(cx, cy).getMaximumBin()*0.04-3),
										(histGroups_X_epg_E.getItem(cx, cy).getMaximumBin()*0.04-1));
					fXe.setParameter(0, histGroups_X_epg_E.getItem(cx, cy).getMax());
					fXe.setParameter(1, (histGroups_X_epg_E.getItem(cx, cy).getMaximumBin()*0.04-2));
					fXe.setParameter(2, histGroups_X_epg_E.getItem(cx, cy).getRMS()/2);
					fXe.setParLimits(2, histGroups_X_epg_E.getItem(cx, cy).getRMS()/4,
										histGroups_X_epg_E.getItem(cx, cy).getRMS());
					DataFitter.fit(fXe, histGroups_X_epg_E.getItem(cx, cy), "Q");
					mean_Xe.add(fXe.getParameter(1), cx, cy);
					sigma_Xe.add(Math.abs(fXe.getParameter(2)), cx, cy);
				}
			}
		}
		
		IndexedList<GraphErrors> htheta_cone_gamma = new IndexedList<GraphErrors>(1);
		IndexedList<GraphErrors> hX_eg_m = new IndexedList<GraphErrors>(1);
		IndexedList<GraphErrors> hX_epg_pt = new IndexedList<GraphErrors>(1);
		IndexedList<GraphErrors> hX_epg_E = new IndexedList<GraphErrors>(1);
		
		for(int cx = 1; cx <= 22; cx++)
		{
			GraphErrors htheta_cone_gamma_ini = new GraphErrors();
			htheta_cone_gamma.add(htheta_cone_gamma_ini, cx);
			GraphErrors hX_eg_m_ini = new GraphErrors();
			hX_eg_m.add(hX_eg_m_ini, cx);
			GraphErrors hX_epg_pt_ini = new GraphErrors();
			hX_epg_pt.add(hX_epg_pt_ini, cx);
			GraphErrors hX_epg_E_ini = new GraphErrors();
			hX_epg_E.add(hX_epg_E_ini, cx);
			for(int cy = 0; cy <=22; cy++)
			{
				if(mean_tc.hasItem(cx, cy) && sigma_tc.hasItem(cx, cy))
				{
					htheta_cone_gamma.getItem(cx).addPoint(cy, mean_tc.getItem(cx, cy), 0, sigma_tc.getItem(cx, cy));
				}
				if(mean_pm.hasItem(cx, cy) && sigma_pm.hasItem(cx, cy))
				{
					hX_eg_m.getItem(cx).addPoint(cy, mean_pm.getItem(cx, cy), 0, sigma_pm.getItem(cx, cy));
				}
				if(mean_Xpt.hasItem(cx, cy) && sigma_Xpt.hasItem(cx, cy))
				{
					hX_epg_pt.getItem(cx).addPoint(cy, mean_Xpt.getItem(cx, cy), 0, sigma_Xpt.getItem(cx, cy));
				}
				if(mean_Xe.hasItem(cx, cy) && sigma_Xe.hasItem(cx, cy))
				{
					hX_epg_E.getItem(cx).addPoint(cy, mean_Xe.getItem(cx, cy), 0, sigma_Xe.getItem(cx, cy));
				}
			}
		}
		
		for(int cx = 1; cx <= 22; cx++)
		{
			JFrame framecxX = new JFrame("FTCAL Crystalized Column x = " + cx + " Exclusivity Variables");
			framecxX.setSize(1000, 1000);
			EmbeddedCanvas cancxX = new EmbeddedCanvas();
			framecxX.add(cancxX);
			framecxX.setLocationRelativeTo(null);
			framecxX.setVisible(true);
			cancxX.divide(2, 2);
			cancxX.cd(0);
			cancxX.setFont("Arial");
			htheta_cone_gamma.getItem(cx).setTitle("#Delta#theta_cone(#gamma)");
			htheta_cone_gamma.getItem(cx).setTitleX("Crystalized Row y");
			htheta_cone_gamma.getItem(cx).setTitleY("#mu +/- #sigma [#degree]");
			htheta_cone_gamma.getItem(cx).setMarkerSize(5);
			htheta_cone_gamma.getItem(cx).setMarkerStyle(1);
			htheta_cone_gamma.getItem(cx).setLineThickness(1);
			cancxX.getPad(0).setTitleFontSize(32);
			cancxX.getPad(0).setAxisTitleFontSize(32);
			cancxX.getPad(0).setAxisLabelFontSize(24);
			cancxX.getPad(0).setAxisRange(0, 24, -3, 3);
			cancxX.getPad(0).setStatBoxFontSize(18);
			cancxX.draw(htheta_cone_gamma.getItem(cx), "same");
			cancxX.cd(1);
			cancxX.setFont("Arial");
			hX_eg_m.getItem(cx).setTitle("M_X_(ep->e'#gamma)");
			hX_eg_m.getItem(cx).setTitleX("Crystalized Row y");
			hX_eg_m.getItem(cx).setTitleY("#mu +/- #sigma [GeV]");
			hX_eg_m.getItem(cx).setMarkerSize(5);
			hX_eg_m.getItem(cx).setMarkerStyle(1);
			hX_eg_m.getItem(cx).setLineThickness(1);
			cancxX.getPad(1).setTitleFontSize(32);
			cancxX.getPad(1).setAxisTitleFontSize(32);
			cancxX.getPad(1).setAxisLabelFontSize(24);
			cancxX.getPad(1).setAxisRange(0, 24, 0, 2);
			cancxX.getPad(1).setStatBoxFontSize(18);
			cancxX.draw(hX_eg_m.getItem(cx), "same");
			cancxX.cd(2);
			cancxX.setFont("Arial");
			hX_epg_pt.getItem(cx).setTitle("#p_#rho_X_(ep->e'p'#gamma)");
			hX_epg_pt.getItem(cx).setTitleX("Crystalized Row y");
			hX_epg_pt.getItem(cx).setTitleY("#mu +/- #sigma [GeV]");
			hX_epg_pt.getItem(cx).setMarkerSize(5);
			hX_epg_pt.getItem(cx).setMarkerStyle(1);
			hX_epg_pt.getItem(cx).setLineThickness(1);
			cancxX.getPad(2).setTitleFontSize(32);
			cancxX.getPad(2).setAxisTitleFontSize(32);
			cancxX.getPad(2).setAxisLabelFontSize(24);
			cancxX.getPad(2).setAxisRange(0, 24, -0.5, 0.5);
			cancxX.getPad(2).setStatBoxFontSize(18);
			cancxX.draw(hX_epg_pt.getItem(cx), "same");
			cancxX.cd(3);
			cancxX.setFont("Arial");
			hX_epg_E.getItem(cx).setTitle("E_X_(ep->e'p'#gamma)");
			hX_epg_E.getItem(cx).setTitleX("Crystalized Row y");
			hX_epg_E.getItem(cx).setTitleY("#mu +/- #sigma [GeV]");
			hX_epg_E.getItem(cx).setMarkerSize(5);
			hX_epg_E.getItem(cx).setMarkerStyle(1);
			hX_epg_E.getItem(cx).setLineThickness(1);
			cancxX.getPad(3).setTitleFontSize(32);
			cancxX.getPad(3).setAxisTitleFontSize(32);
			cancxX.getPad(3).setAxisLabelFontSize(24);
			cancxX.getPad(3).setAxisRange(0, 24, -1, 1);
			cancxX.getPad(3).setStatBoxFontSize(18);
			cancxX.draw(hX_epg_E.getItem(cx), "same");
		}
		
	//	Set crystallized coordinates (cx[cb], cy[cb])
		
		int[] cx = new int[] {7,7,6,4,5};
		int[] cy = new int[] {15,16,15,12};
		
		for(int cb = 0; cb < cx.length ; cb++)
		{
			JFrame frameX = new JFrame("FTCAL Crystalized Coordinate (x,y) = (" + cx[cb] + "," + cy[cb] +") Exclusivity Variables");
			frameX.setSize(1000, 1000);
			EmbeddedCanvas canX = new EmbeddedCanvas();
			frameX.add(canX);
			frameX.setLocationRelativeTo(null);
			frameX.setVisible(true);
			canX.divide(2, 2);
			canX.cd(0);
			canX.setFont("Arial");
			histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).setTitle("#Delta#theta_cone(#gamma)");
			histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).setTitleX("#Delta#theta_cone(#gamma) [#degree]");
			histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).setTitleY("Counts");
			histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).setOptStat(10);
			canX.getPad(0).setTitleFontSize(32);
			canX.getPad(0).setAxisTitleFontSize(32);
			canX.getPad(0).setAxisLabelFontSize(24);
			canX.getPad(0).setStatBoxFontSize(18);
			canX.draw(histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]), "same");
			F1D ftc = new F1D("ftc", "[amp]*gaus(x,[mean],[sigma])",
								histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).getMaximumBin()*0.05, 
								histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).getRMS()*3);
			ftc.setParameter(0, histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).getMax());
			ftc.setParLimits(0, histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).getMax()+1,
								histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).getMax()-1);
			ftc.setParameter(1, (histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).getMaximumBin()*0.05));
			ftc.setParLimits(1, ((histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).getMaximumBin()-1)*0.05),
								((histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).getMaximumBin()+1)*0.05));
			ftc.setParameter(2, histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).getRMS()/2);
			ftc.setParLimits(2, histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).getRMS()/4,
								histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]).getRMS());
			DataFitter.fit(ftc, histGroups_theta_cone_gamma.getItem(cx[cb], cy[cb]), "Q");
			ftc.setLineColor(2);
			ftc.setLineWidth(3);
			ftc.setOptStat(1100);
			canX.draw(ftc, "same");
			canX.cd(1);
			canX.setFont("Arial");
			histGroups_X_eg_m.getItem(cx[cb], cy[cb]).setTitle("M_X_(ep->e'#gamma)");
			histGroups_X_eg_m.getItem(cx[cb], cy[cb]).setTitleX("M_X [GeV]");
			histGroups_X_eg_m.getItem(cx[cb], cy[cb]).setTitleY("Counts");
			histGroups_X_eg_m.getItem(cx[cb], cy[cb]).setOptStat(10);
			canX.getPad(1).setTitleFontSize(32);
			canX.getPad(1).setAxisTitleFontSize(32);
			canX.getPad(1).setAxisLabelFontSize(24);
			canX.getPad(1).setStatBoxFontSize(18);
			canX.draw(histGroups_X_eg_m.getItem(cx[cb], cy[cb]), "same");
			F1D fpm = new F1D("fpm", "[amp]*gaus(x,[mean],[sigma])", (histGroups_X_eg_m.getItem(cx[cb], cy[cb]).getMaximumBin()*0.03-1),
								(histGroups_X_eg_m.getItem(cx[cb], cy[cb]).getMaximumBin()*0.03+1));
			fpm.setParameter(0, histGroups_X_eg_m.getItem(cx[cb], cy[cb]).getMax());
			fpm.setParameter(1, histGroups_X_eg_m.getItem(cx[cb], cy[cb]).getMaximumBin()*0.03);
			fpm.setParameter(2, histGroups_X_eg_m.getItem(cx[cb], cy[cb]).getRMS()/2);
			fpm.setParLimits(2, histGroups_X_eg_m.getItem(cx[cb], cy[cb]).getRMS()/4,
								histGroups_X_eg_m.getItem(cx[cb], cy[cb]).getRMS());
			DataFitter.fit(fpm, histGroups_X_eg_m.getItem(cx[cb], cy[cb]), "Q");
			fpm.setLineColor(2);
			fpm.setLineWidth(3);
			fpm.setOptStat(1100);
			canX.draw(fpm, "same");
			canX.cd(2);
			canX.setFont("Arial");
			histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).setTitle("p_#rho_X_(ep->e'p'#gamma)");
			histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).setTitleX("p_#rho_X [GeV]");
			histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).setTitleY("Counts");
			histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).setOptStat(10);
			canX.getPad(2).setTitleFontSize(32);
			canX.getPad(2).setAxisTitleFontSize(32);
			canX.getPad(2).setAxisLabelFontSize(24);
			canX.getPad(2).setStatBoxFontSize(18);
			canX.draw(histGroups_X_epg_pt.getItem(cx[cb], cy[cb]), "same");
			F1D fXpt = new F1D("fXpt", "[amp]*gaus(x,[mean],[sigma])",
								histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).getMaximumBin()*0.008,
								histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).getRMS()*3);
			fXpt.setParameter(0, histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).getMax());
			fXpt.setParLimits(0, histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).getMax()+1,
								histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).getMax()-1);
			fXpt.setParameter(1, (histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).getMaximumBin()*0.008));
			fXpt.setParLimits(1, ((histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).getMaximumBin()-1)*0.008),
								((histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).getMaximumBin()+1)*0.008));
			fXpt.setParameter(2, histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).getRMS()/2);
			fXpt.setParLimits(2, histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).getRMS()/4,
								histGroups_X_epg_pt.getItem(cx[cb], cy[cb]).getRMS());
			DataFitter.fit(fXpt, histGroups_X_epg_pt.getItem(cx[cb], cy[cb]), "Q");
			fXpt.setLineColor(2);
			fXpt.setLineWidth(3);
			fXpt.setOptStat(1100);
			canX.draw(fXpt, "same");
			canX.cd(3);
			canX.setFont("Arial");
			histGroups_X_epg_E.getItem(cx[cb], cy[cb]).setTitle("E_X_(ep->e'p'#gamma)");
			histGroups_X_epg_E.getItem(cx[cb], cy[cb]).setTitleX("E_X [GeV]");
			histGroups_X_epg_E.getItem(cx[cb], cy[cb]).setTitleY("Counts");
			histGroups_X_epg_E.getItem(cx[cb], cy[cb]).setOptStat(10);
			canX.getPad(3).setTitleFontSize(32);
			canX.getPad(3).setAxisTitleFontSize(32);
			canX.getPad(3).setAxisLabelFontSize(24);
			canX.getPad(3).setStatBoxFontSize(18);
			canX.draw(histGroups_X_epg_E.getItem(cx[cb], cy[cb]), "same");
			F1D fXe = new F1D("fXe", "[amp]*gaus(x,[mean],[sigma])", (histGroups_X_epg_E.getItem(cx[cb], cy[cb]).getMaximumBin()*0.04-3),
								(histGroups_X_epg_E.getItem(cx[cb], cy[cb]).getMaximumBin()*0.04-1));
			fXe.setParameter(0, histGroups_X_epg_E.getItem(cx[cb], cy[cb]).getMax());
			fXe.setParameter(1, (histGroups_X_epg_E.getItem(cx[cb], cy[cb]).getMaximumBin()*0.04-2));
			fXe.setParameter(2, histGroups_X_epg_E.getItem(cx[cb], cy[cb]).getRMS()/2);
			fXe.setParLimits(2, histGroups_X_epg_E.getItem(cx[cb], cy[cb]).getRMS()/4,
								histGroups_X_epg_E.getItem(cx[cb], cy[cb]).getRMS());
			DataFitter.fit(fXe, histGroups_X_epg_E.getItem(cx[cb], cy[cb]), "Q");
			fXe.setLineColor(2);
			fXe.setLineWidth(3);
			fXe.setOptStat(1100);
			canX.draw(fXe, "same");
		}
	}
}