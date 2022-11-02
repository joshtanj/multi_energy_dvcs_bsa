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

public class dc_theta_phi_binner {
	
	static HipoDataSource reader = new HipoDataSource();
	static H2F hdc_det_y_vs_x = new H2F("dc_det_y_vs_x", "dc_det_y_vs_x", 500, -125, 125, 500, -125, 125);
	static H1F hdc_det_theta = new H1F("dc_det_theta", "dc_det_theta", 200, 3, 33);
	static H2F hdc_det_theta_vs_phi = new H2F("dc_det_theta_vs_phi", "dc_det_theta_vs_phi", 1500, -150, 210, 500, 3, 33);
	static H2F hdc_det_theta_vs_phi_rot = new H2F("dc_det_theta_vs_phi_rot", "dc_det_theta_vs_phi_rot", 500, -60, 60, 500, 3, 33);
	static H2F hdc_det_theta_vs_phi_rot_mir = new H2F("dc_det_theta_vs_phi_rot_mir", "dc_det_theta_vs_phi_rot_mir", 500, -60, 60, 500, 3, 33);
	static H1F hdc_det_phi_rot_mir_prj = new H1F("dc_det_phi_rot_mir_prj", "dc_det_phi_rot_mir_prj", 200, -30, 30);
	
	static IndexedList<H1F> histGroups_dc_det_theta_bin= new IndexedList<H1F>(1);
	public static void dc_det_theta_bin_histos() {
		for(int ibin = 0; ibin < 9; ibin++) {
			H1F dc_det_theta_bin = new H1F("dc_det_theta_bin", "dc_det_theta_bin", 240, 6, 33);
			histGroups_dc_det_theta_bin.add(dc_det_theta_bin, ibin);
		}
	}
	
	static H1F hdc_det_phi_rot_mir_prj_bin1 = new H1F("dc_det_phi_rot_mir_prj_bin1", "dc_det_phi_rot_mir_prj_bin1", 200, -30, 30);
	static H1F hdc_det_phi_rot_mir_prj_bin2 = new H1F("dc_det_phi_rot_mir_prj_bin2", "dc_det_phi_rot_mir_prj_bin2", 200, -30, 30);
	
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
			float theta_deg_e = 0;
			float phi_deg_e = -200;
			float phi_rot_deg_e = -200;
			double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0.000511);
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0.938272);
		//	double E = 6.535;
			double E = 7.546;
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
				hdc_det_y_vs_x.fill(dc_hit_e.x(), dc_hit_e.y());
				hdc_det_theta.fill(theta_deg_e);
				if(theta_deg_e > 7 && theta_deg_e < 33)
				{
				hdc_det_theta_vs_phi.fill(phi_deg_e, theta_deg_e);
				hdc_det_theta_vs_phi_rot.fill(phi_rot_deg_e, theta_deg_e);
				hdc_det_theta_vs_phi_rot_mir.fill(Math.abs(phi_rot_deg_e), theta_deg_e);
				hdc_det_phi_rot_mir_prj.fill(Math.abs(phi_rot_deg_e));
				}
				double[] theta_bnd  = new double[] {7, 8.5, 10, 12, 14, 16, 18, 21, 24, 33};
				for(int thbin = 0; thbin < 9; thbin++)
				{
					if(theta_deg_e > theta_bnd[thbin] && theta_deg_e <= theta_bnd[thbin+1])
					{
						histGroups_dc_det_theta_bin.getItem(thbin).fill(theta_deg_e);
						break;
					}
				}
			}
		}
	}
	
	public static void main(String[] args) {
		
		dc_det_theta_bin_histos();
	
	//	reader.open("C:/Users/joshtanj/Documents/download/merged_skim_e_bank_6535MeV_skim_elastic_5950_5899_5963_5958_5885.hipo");
		reader.open("C:/Users/joshtanj/Documents/download/merged_skim_ep_bank_7546MeV_skim_elastic.hipo");
		
		int eventCounter = 1;
		while(reader.hasEvent())// && eventCounter < 1000000)
		{
			processEvent(reader.getNextEvent(), eventCounter);
			if(eventCounter%50000 == 0) System.out.println("Event: " + eventCounter);
			eventCounter++;
		}
		int Nevent = eventCounter-1;
		System.out.println("Number of events: " + Nevent);
		
		JFrame frame_dc_hit = new JFrame("DC Hits");
		frame_dc_hit.setSize(500, 500);
		EmbeddedCanvas can_dc_hit = new EmbeddedCanvas();
		frame_dc_hit.add(can_dc_hit);
		frame_dc_hit.setLocationRelativeTo(null);
		frame_dc_hit.setVisible(true);
		can_dc_hit.divide(1, 1);
		can_dc_hit.cd(0);
		can_dc_hit.setFont("Arial");
		hdc_det_y_vs_x.setTitle("DC y vs. x");
		hdc_det_y_vs_x.setTitleX("x [cm]");
		hdc_det_y_vs_x.setTitleY("y [cm]");
		can_dc_hit.getPad(0).setTitleFontSize(32);
		can_dc_hit.getPad(0).setAxisTitleFontSize(32);
		can_dc_hit.getPad(0).setAxisLabelFontSize(24);
		can_dc_hit.getPad(0).setStatBoxFontSize(18);
		can_dc_hit.draw(hdc_det_y_vs_x, "same");
		
		JFrame frame_dc_hit_theta = new JFrame("DC Hits #theta");
		frame_dc_hit_theta.setSize(500, 500);
		EmbeddedCanvas can_dc_hit_theta = new EmbeddedCanvas();
		frame_dc_hit_theta.add(can_dc_hit_theta);
		frame_dc_hit_theta.setLocationRelativeTo(null);
		frame_dc_hit_theta.setVisible(true);
		can_dc_hit_theta.divide(1, 1);
		can_dc_hit_theta.cd(0);
		can_dc_hit_theta.setFont("Arial");
		hdc_det_theta.setTitle("DC #theta");
		hdc_det_theta.setTitleX("#theta [#degree]");
		hdc_det_theta.setTitleY("Counts");
		hdc_det_theta.setOptStat(10);
		can_dc_hit_theta.getPad(0).setTitleFontSize(32);
		can_dc_hit_theta.getPad(0).setAxisTitleFontSize(32);
		can_dc_hit_theta.getPad(0).setAxisLabelFontSize(24);
		can_dc_hit_theta.getPad(0).setStatBoxFontSize(18);
		can_dc_hit_theta.draw(hdc_det_theta, "same");
		
		JFrame frame_dc_hit_theta_vs_phi = new JFrame("DC Hits #theta vs. #phi");
		frame_dc_hit_theta_vs_phi.setSize(1500, 500);
		EmbeddedCanvas can_dc_hit_theta_vs_phi = new EmbeddedCanvas();
		frame_dc_hit_theta_vs_phi.add(can_dc_hit_theta_vs_phi);
		frame_dc_hit_theta_vs_phi.setLocationRelativeTo(null);
		frame_dc_hit_theta_vs_phi.setVisible(true);
		can_dc_hit_theta_vs_phi.divide(1, 1);
		can_dc_hit_theta_vs_phi.cd(0);
		can_dc_hit_theta_vs_phi.setFont("Arial");
		hdc_det_theta_vs_phi.setTitle("DC #theta vs #phi");
		hdc_det_theta_vs_phi.setTitleX("#phi [#degree]");
		hdc_det_theta_vs_phi.setTitleY("#theta [#degree]");
		can_dc_hit_theta_vs_phi.getPad(0).setTitleFontSize(32);
		can_dc_hit_theta_vs_phi.getPad(0).setAxisTitleFontSize(32);
		can_dc_hit_theta_vs_phi.getPad(0).setAxisLabelFontSize(24);
		can_dc_hit_theta_vs_phi.getPad(0).setStatBoxFontSize(18);
		can_dc_hit_theta_vs_phi.draw(hdc_det_theta_vs_phi, "same");
		
		JFrame frame_dc_hit_theta_vs_phi_prj = new JFrame("DC Hits #phi Projections");
		frame_dc_hit_theta_vs_phi_prj.setSize(1500, 500);
		EmbeddedCanvas can_dc_hit_theta_vs_phi_prj = new EmbeddedCanvas();
		frame_dc_hit_theta_vs_phi_prj.add(can_dc_hit_theta_vs_phi_prj);
		frame_dc_hit_theta_vs_phi_prj.setLocationRelativeTo(null);
		frame_dc_hit_theta_vs_phi_prj.setVisible(true);
		can_dc_hit_theta_vs_phi_prj.divide(3, 1);
		can_dc_hit_theta_vs_phi_prj.cd(0);
		can_dc_hit_theta_vs_phi_prj.setFont("Arial");
		hdc_det_theta_vs_phi_rot.setTitle("Superpositioned DC #theta vs. #phi");
		hdc_det_theta_vs_phi_rot.setTitleX("#phi [#degree]");
		hdc_det_theta_vs_phi_rot.setTitleY("#theta [#degree]");
		can_dc_hit_theta_vs_phi_prj.getPad(0).setTitleFontSize(32);
		can_dc_hit_theta_vs_phi_prj.getPad(0).setAxisTitleFontSize(32);
		can_dc_hit_theta_vs_phi_prj.getPad(0).setAxisLabelFontSize(24);
		can_dc_hit_theta_vs_phi_prj.getPad(0).setStatBoxFontSize(18);
		can_dc_hit_theta_vs_phi_prj.draw(hdc_det_theta_vs_phi_rot, "same");
		can_dc_hit_theta_vs_phi_prj.cd(1);
		can_dc_hit_theta_vs_phi_prj.setFont("Arial");
		hdc_det_theta_vs_phi_rot_mir.setTitle("Superpositioned Mirrored DC #theta vs. #phi");
		hdc_det_theta_vs_phi_rot_mir.setTitleX("#phi [#degree]");
		hdc_det_theta_vs_phi_rot_mir.setTitleY("#theta [#degree]");
		can_dc_hit_theta_vs_phi_prj.getPad(1).setTitleFontSize(32);
		can_dc_hit_theta_vs_phi_prj.getPad(1).setAxisTitleFontSize(32);
		can_dc_hit_theta_vs_phi_prj.getPad(1).setAxisLabelFontSize(24);
		can_dc_hit_theta_vs_phi_prj.getPad(1).setStatBoxFontSize(18);
		can_dc_hit_theta_vs_phi_prj.draw(hdc_det_theta_vs_phi_rot_mir, "same");
		can_dc_hit_theta_vs_phi_prj.cd(2);
		can_dc_hit_theta_vs_phi_prj.setFont("Arial");
		hdc_det_phi_rot_mir_prj.setTitle("Superpositioned Mirrored DC #phi");
		hdc_det_phi_rot_mir_prj.setTitleX("#phi [#degree]");
		hdc_det_phi_rot_mir_prj.setTitleY("Counts");
		hdc_det_phi_rot_mir_prj.setOptStat(10);
		can_dc_hit_theta_vs_phi_prj.getPad(2).setTitleFontSize(32);
		can_dc_hit_theta_vs_phi_prj.getPad(2).setAxisTitleFontSize(32);
		can_dc_hit_theta_vs_phi_prj.getPad(2).setAxisLabelFontSize(24);
		can_dc_hit_theta_vs_phi_prj.getPad(2).setStatBoxFontSize(18);
		can_dc_hit_theta_vs_phi_prj.draw(hdc_det_phi_rot_mir_prj, "same");
		
		double[] theta_bnd  = new double[] {7, 8.5, 10, 12, 14, 16, 18, 21, 24, 33};
		
		JFrame frame_dc_hit_theta_bin = new JFrame("DC Hits #theta Bins");
		frame_dc_hit_theta_bin.setSize(1500, 1000);
		EmbeddedCanvas can_dc_hit_theta_bin = new EmbeddedCanvas();
		frame_dc_hit_theta_bin.add(can_dc_hit_theta_bin);
		frame_dc_hit_theta_bin.setLocationRelativeTo(null);
		frame_dc_hit_theta_bin.setVisible(true);
		can_dc_hit_theta_bin.divide(3, 3);	
		for(int thbin = 0; thbin < 9; thbin++)
		{
			can_dc_hit_theta_bin.cd(thbin);
			can_dc_hit_theta_bin.setFont("Arial");
			histGroups_dc_det_theta_bin.getItem(thbin).setTitle("DC " + theta_bnd[thbin] + "#degree < #theta <= "
																	+ theta_bnd[thbin+1] + "#degree");
			histGroups_dc_det_theta_bin.getItem(thbin).setTitleX("#theta [#degree]");
			histGroups_dc_det_theta_bin.getItem(thbin).setTitleY("Counts");
			histGroups_dc_det_theta_bin.getItem(thbin).setOptStat(110);
			can_dc_hit_theta_bin.getPad(thbin).setTitleFontSize(32);
			can_dc_hit_theta_bin.getPad(thbin).setAxisTitleFontSize(32);
			can_dc_hit_theta_bin.getPad(thbin).setAxisLabelFontSize(24);
			can_dc_hit_theta_bin.getPad(thbin).setStatBoxFontSize(18);
			can_dc_hit_theta_bin.draw(histGroups_dc_det_theta_bin.getItem(thbin), "same");
		}
		
		for(int bini = 0; bini < 200; bini++){
			if(hdc_det_phi_rot_mir_prj_bin1.getEntries() < hdc_det_phi_rot_mir_prj.getEntries()/2){
				for(int filli = 0; filli < hdc_det_phi_rot_mir_prj.getBinContent(bini); filli++){
					hdc_det_phi_rot_mir_prj_bin1.fill(-30+0.3*bini);
				}
			}
			else{
				for(int filli = 0; filli < hdc_det_phi_rot_mir_prj.getBinContent(bini); filli++){
					hdc_det_phi_rot_mir_prj_bin2.fill(-30+0.3*bini);
				}
			}
		}
		
		JFrame frame_dc_hit_phi_bin = new JFrame("DC Hits #phi Bins");
		frame_dc_hit_phi_bin.setSize(1000, 500);
		EmbeddedCanvas can_dc_hit_phi_bin = new EmbeddedCanvas();
		frame_dc_hit_phi_bin.add(can_dc_hit_phi_bin);
		frame_dc_hit_phi_bin.setLocationRelativeTo(null);
		frame_dc_hit_phi_bin.setVisible(true);
		can_dc_hit_phi_bin.divide(2, 1);
		can_dc_hit_phi_bin.cd(0);
		can_dc_hit_phi_bin.setFont("Arial");
		hdc_det_phi_rot_mir_prj_bin1.setTitle("DC #phi Bin 1");
		hdc_det_phi_rot_mir_prj_bin1.setTitleX("#phi [#degree]");
		hdc_det_phi_rot_mir_prj_bin1.setTitleY("Counts");
		hdc_det_phi_rot_mir_prj_bin1.setOptStat(110);
		can_dc_hit_phi_bin.getPad(0).setTitleFontSize(32);
		can_dc_hit_phi_bin.getPad(0).setAxisTitleFontSize(32);
		can_dc_hit_phi_bin.getPad(0).setAxisLabelFontSize(24);
		can_dc_hit_phi_bin.getPad(0).setStatBoxFontSize(18);
		can_dc_hit_phi_bin.draw(hdc_det_phi_rot_mir_prj_bin1, "same");
		can_dc_hit_phi_bin.cd(1);
		can_dc_hit_phi_bin.setFont("Arial");
		hdc_det_phi_rot_mir_prj_bin2.setTitle("DC #phi Bin 2");
		hdc_det_phi_rot_mir_prj_bin2.setTitleX("#phi [#degree]");
		hdc_det_phi_rot_mir_prj_bin2.setTitleY("Counts");
		hdc_det_phi_rot_mir_prj_bin2.setOptStat(110);
		can_dc_hit_phi_bin.getPad(1).setTitleFontSize(32);
		can_dc_hit_phi_bin.getPad(1).setAxisTitleFontSize(32);
		can_dc_hit_phi_bin.getPad(1).setAxisLabelFontSize(24);
		can_dc_hit_phi_bin.getPad(1).setStatBoxFontSize(18);
		can_dc_hit_phi_bin.draw(hdc_det_phi_rot_mir_prj_bin2, "same");
	}
}