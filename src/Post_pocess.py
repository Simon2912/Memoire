import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import math

class StructureVisualizer:
    def __init__(self, coords, U, connectivity, elem_types, initial_point_index, convergence, p_bsc):
        self.coords = coords
        self.U = U
        self.connectivity = connectivity
        self.elem_types = elem_types
        self.initial_point_index = initial_point_index
        self.convergence = convergence
        self.p_bsc = p_bsc
        
        # --- CORRECTION 1 : Arrondi propre avec NumPy ---
        if self.p_bsc is not None:
            self.p_bsc = np.round(self.p_bsc, 1)
        # Détection automatique du nombre de DDL
        n_nodes = len(coords[0])
        n_dof_total = len(U)
        dof_per_node = n_dof_total // n_nodes 
        
        self.Ux = U[0::dof_per_node]
        self.Uy = U[1::dof_per_node]
        self.Uz = U[2::dof_per_node]
        
        self.x0 = coords[0]
        self.y0 = coords[1]
        self.z0 = coords[2]
        
        self.show_nodes = False
        self.selected_effort = 'Aucun' 
        self.scale = 1.0
        self.color_map = {1: 'green', 2: 'red', 3: 'blue'}
        
        self.fig = plt.figure(figsize=(12, 7)) 
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        # On ajuste l'espace principal pour laisser de la place à droite (jusqu'à 72% de la largeur)
        plt.subplots_adjust(bottom=0.25, right=0.72) 
        
        # --- CORRECTION DU BUG DE RÉTRÉCISSEMENT ---
        # On crée un axe fixe dédié UNIQUEMENT à la Colorbar.
        # Format : [gauche, bas, largeur, hauteur]
        self.cax = self.fig.add_axes([0.74, 0.25, 0.02, 0.5]) 
        self.cax.set_visible(False) # Caché par défaut
        self.cbar = None
        # -------------------------------------------
        
        self._init_widgets()
        self.print_stats()
        self.print_efforts_table() 
        self.update_plot()

    def print_stats(self):
        print(f"\n=== STATISTIQUES DE DÉPLACEMENT ===")
        print(f"Deplacement max en x : {np.max(np.abs(self.Ux)):.4f} m")
        print(f"Deplacement max en y : {np.max(np.abs(self.Uy)):.4f} m")
        print(f"Deplacement max en z : {np.max(np.abs(self.Uz)):.10f} m")
        
        indices = np.atleast_1d(self.initial_point_index)
        for i, idx in enumerate(indices):
            print(f"Deplacement du noeud initial (noeud {idx+1}): "
                  f"Ux={self.Ux[idx]:.4f}, Uy={self.Uy[idx]:.4f}, Uz={self.Uz[idx]:.4f}")

    def print_efforts_table(self):
        if self.p_bsc is None:
            print("\nERREUR : 'p_bsc' manquant. Impossible d'afficher les efforts.")
            return
            
        print("\n=== TABLEAU DES EFFORTS INTERNES ===")
        for idx in range(len(self.connectivity)):
            N = self.p_bsc[idx][0]
            My = max(abs(self.p_bsc[idx][1]), abs(self.p_bsc[idx][2]))
            Mz = max(abs(self.p_bsc[idx][3]), abs(self.p_bsc[idx][4]))
            M = math.sqrt(My**2 + Mz**2)
            Mt = self.p_bsc[idx][5]
            
            print(f"Élém {idx:2d} | Axial(N): {N:10.2f} N | Flex(My): {My:10.2f} Nm | Flex(Mz): {Mz:10.2f} Nm | Flex(M): {M:10.2f} Nm | Torsion(Mt): {Mt:10.2f} Nm")
        print("====================================\n")
    
    def _init_widgets(self):
        ax_slider = plt.axes([0.25, 0.1, 0.45, 0.03])
        self.slider = Slider(ax_slider, 'Échelle', 0.0, 5.0, valinit=self.scale)
        self.slider.on_changed(self._on_slider_change)

        ax_btn_nodes = plt.axes([0.82, 0.04, 0.15, 0.04])
        self.btn_nodes = Button(ax_btn_nodes, 'N° Nœuds', color='lightgray', hovercolor='0.85')
        self.btn_nodes.on_clicked(self._toggle_labels)

        ax_radio = plt.axes([0.82, 0.12, 0.15, 0.22], facecolor='lightgray')
        self.radio_efforts = RadioButtons(ax_radio, ('Aucun', 'Normal (N)', 'Flexion (My)', 'Flexion (Mz)', 'Flexion (M)', 'Torsion (Mt)'))
        self.radio_efforts.on_clicked(self._on_effort_select)

        ax_btn_reset = plt.axes([0.25, 0.04, 0.12, 0.04])
        self.btn_reset = Button(ax_btn_reset, 'Reset', color='lightgray', hovercolor='0.85')
        self.btn_reset.on_clicked(self._reset_view)

    def _on_slider_change(self, val):
        self.scale = val
        self.update_plot()

    def _toggle_labels(self, event):
        self.show_nodes = not self.show_nodes
        self.update_plot()

    def _on_effort_select(self, label):
        self.selected_effort = label
        self.update_plot()

    def _reset_view(self, event):
        self.slider.reset()
        
    def update_plot(self):
        self.ax.clear()
        
        xd = self.x0 + self.scale * self.Ux
        yd = self.y0 + self.scale * self.Uy
        zd = self.z0 + self.scale * self.Uz
        
        cmap = cm.jet
        norm = None
        efforts_list = []

        # --- GESTION DE LA COLORBAR FIXE ---
        if self.selected_effort != 'Aucun' and self.p_bsc is not None:
            for idx in range(len(self.connectivity)):
                # Extraction de la valeur selon le choix
                if self.selected_effort == 'Normal (N)':
                    val = self.p_bsc[idx][0]
                elif self.selected_effort == 'Flexion (My)':
                    val = max(abs(self.p_bsc[idx][1]), abs(self.p_bsc[idx][2]))
                elif self.selected_effort == 'Flexion (Mz)':
                    val = max(abs(self.p_bsc[idx][3]), abs(self.p_bsc[idx][4]))
                elif self.selected_effort == 'Flexion (M)':
                    val = math.sqrt(max(abs(self.p_bsc[idx][1]), abs(self.p_bsc[idx][2]))**2 + max(abs(self.p_bsc[idx][3]), abs(self.p_bsc[idx][4]))**2)
                elif self.selected_effort == 'Torsion (Mt)':
                    val = self.p_bsc[idx][5]
                
                efforts_list.append(val)
                # ARRONDISSEMENT À 2 DÉCIMALES (Suppression du bruit numérique)

            if len(efforts_list) > 0:
                vmin, vmax = min(efforts_list), max(efforts_list)
                
                # --- ANCRAGE DE L'ÉCHELLE À ZÉRO ---
                if vmin > 0:
                    vmin = 0.0  # Que du positif : on force le minimum à 0
                elif vmax < 0:
                    vmax = 0.0  # Que du négatif : on force le maximum à 0
                elif vmin == 0 and vmax == 0:
                    vmin = -1.0; vmax = 1.0 # Cas particulier où absolument tout est nul
                    
                # --- CHOIX DE LA NORME ET DE LA PALETTE ---
                if vmin < 0 and vmax > 0:
                    norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=vmin, vmax=vmax)
                    cmap = cm.jet
                else:
                    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
                    cmap = cm.jet 
                    
                sm = cm.ScalarMappable(norm=norm, cmap=cmap)
                sm.set_array([])
                
                self.cax.clear()
                self.cbar = self.fig.colorbar(sm, cax=self.cax)
                
                # --- DÉSACTIVATION DE L'OFFSET ---
                self.cbar.formatter.set_useOffset(False) 
                self.cbar.formatter.set_powerlimits((-3, 4))
                self.cbar.update_ticks()
                
                self.cbar.set_label(f"Intensité - {self.selected_effort}")
                self.cax.set_visible(True)
            else:
                self.cax.set_visible(False)
        else:
            self.cax.set_visible(False)

        # --- TRACÉ DES ÉLÉMENTS ---
        for idx, (i, j) in enumerate(self.connectivity):
            self.ax.plot(
                [self.x0[i], self.x0[j]], [self.y0[i], self.y0[j]], [self.z0[i], self.z0[j]],
                linestyle='--', color='gray', linewidth=0.5, alpha=0.5
            )
            
            # Application de la couleur en fonction de l'effort
            if self.selected_effort != 'Aucun' and norm is not None and len(efforts_list) > 0:
                color = cmap(norm(efforts_list[idx])) 
                line_width = 3.5 
            else:
                elem_type = self.elem_types[idx]
                color = self.color_map.get(elem_type, 'black') 
                line_width = 2.0

            self.ax.plot(
                [xd[i], xd[j]], [yd[i], yd[j]], [zd[i], zd[j]],
                color=color, linewidth=line_width
            )
        
        # Affichage conditionnel des numéros de nœuds (et leurs déplacements)
        if self.show_nodes:
            target_indices = np.atleast_1d(self.initial_point_index).tolist()
            for n in range(len(self.x0)):
                label_text = f'{n+1}'
                text_color = 'blue'
                if n in target_indices:
                    dx = self.Ux[n]; dy = self.Uy[n]; dz = self.Uz[n]
                    label_text += f"\nUx:{dx:.2e}\nUy:{dy:.2e}\nUz:{dz:.2e}"
                    text_color = 'red'
                self.ax.text(xd[n], yd[n], zd[n], label_text, color=text_color, fontsize=8)

        # Maintien de l'aspect ratio
        all_x = np.concatenate([self.x0, xd])
        all_y = np.concatenate([self.y0, yd])
        all_z = np.concatenate([self.z0, zd])
        
        max_range = np.array([all_x.max()-all_x.min(), all_y.max()-all_y.min(), all_z.max()-all_z.min()]).max() / 2.0
        mid_x = (all_x.max()+all_x.min()) * 0.5
        mid_y = (all_y.max()+all_y.min()) * 0.5
        mid_z = (all_z.max()+all_z.min()) * 0.5

        self.ax.set_xlim(mid_x - max_range, mid_x + max_range)
        self.ax.set_ylim(mid_y - max_range, mid_y + max_range)
        self.ax.set_zlim(mid_z - max_range, mid_z + max_range)

        self.ax.set_xlabel('X [m]')
        self.ax.set_ylabel('Y [m]')
        self.ax.set_zlabel('Z [m]')
        self.ax.set_title("Visualisation de la déformation et des efforts")
      
        info_text = "Convergence" if self.convergence else "Non-convergence"
        self.ax.text2D(0.05, 0.95, info_text, transform=self.ax.transAxes, color='black')
        
        self.fig.canvas.draw_idle()

def post_process(Coord, U, Connect, Elem_Types, initial_point_index, _lambda, p_bsc=None):
    viz = StructureVisualizer(Coord, U, Connect, Elem_Types, initial_point_index, _lambda, p_bsc)
    plt.show()