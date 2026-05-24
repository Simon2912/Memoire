import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from mpl_toolkits.mplot3d import Axes3D

class StructureVisualizer:
    def __init__(self, coords, U, connectivity, elem_types, initial_point_index, convergence):
        self.coords = coords
        self.U = U
        self.connectivity = connectivity
        self.elem_types = elem_types
        self.initial_point_index = initial_point_index
        self.convergence = convergence
        
        # --- CORRECTION ICI ---
        # Détection automatique du nombre de DDL par nœud pour éviter les erreurs
        n_nodes = len(coords[0])
        n_dof_total = len(U)
        dof_per_node = n_dof_total // n_nodes # Devrait donner 6 dans votre cas
        
        # Extraction des déplacements selon le nombre de DDL
        # On prend ux, uy, uz (les 3 premiers) et on saute 'dof_per_node' indices
        self.Ux = U[0::dof_per_node]
        self.Uy = U[1::dof_per_node]
        self.Uz = U[2::dof_per_node]
        # ----------------------
        
        # Coordonnées initiales
        self.x0 = coords[0]
        self.y0 = coords[1]
        self.z0 = coords[2]
        
        # ... (le reste du __init__ reste identique)
        self.show_nodes = False
        self.scale = 1.0
        self.color_map = {1: 'green', 2: 'red', 3: 'blue'}
        self.fig = plt.figure(figsize=(10, 6))
        self.ax = self.fig.add_subplot(111, projection='3d')
        plt.subplots_adjust(bottom=0.25)
        self._init_widgets()
        self.print_stats()
        self.update_plot()

    def print_stats(self):
        """Affiche les statistiques de déplacement dans la console."""
        print(f"Deplacement max en x : {np.max(np.abs(self.Ux)):.4f} m")
        print(f"Deplacement max en y : {np.max(np.abs(self.Uy)):.4f} m")
        print(f"Deplacement max en z : {np.max(np.abs(self.Uz)):.4f} m")

        # Note: self.initial_point_index contient des entiers, on peut itérer directement
        # Assurons-nous que c'est une liste ou un tableau pour éviter des erreurs
        indices = np.atleast_1d(self.initial_point_index)
        for i, idx in enumerate(indices):
            print(f"Deplacement du noeud initial (noeud {idx+1}): "
                  f"Ux={self.Ux[idx]:.4f}, Uy={self.Uy[idx]:.4f}, Uz={self.Uz[idx]:.4f}")

    def _init_widgets(self):
        """Crée les sliders et boutons."""
        # Slider Scale
        ax_slider = plt.axes([0.25, 0.1, 0.45, 0.03])
        self.slider = Slider(ax_slider, 'Échelle', 0.0, 5.0, valinit=self.scale)
        self.slider.on_changed(self._on_slider_change)

        # Bouton Nœuds
        ax_btn_nodes = plt.axes([0.82, 0.1, 0.12, 0.04])
        self.btn_nodes = Button(ax_btn_nodes, 'N° Nœuds', color='lightgray', hovercolor='0.85')
        self.btn_nodes.on_clicked(self._toggle_labels)

        # Bouton Reset
        ax_btn_reset = plt.axes([0.25, 0.05, 0.12, 0.04])
        self.btn_reset = Button(ax_btn_reset, 'Reset', color='lightgray', hovercolor='0.85')
        self.btn_reset.on_clicked(self._reset_view)

    def _on_slider_change(self, val):
        self.scale = val
        self.update_plot()

    def _toggle_labels(self, event):
        self.show_nodes = not self.show_nodes
        self.update_plot()

    def _reset_view(self, event):
        self.slider.reset() # Remet à valinit (1.0)
        
    def update_plot(self):
        """Redessine la structure."""
        self.ax.clear()
        
        # Calcul des coordonnées déformées
        xd = self.x0 + self.scale * self.Ux
        yd = self.y0 + self.scale * self.Uy
        zd = self.z0 + self.scale * self.Uz
        
        # Tracé des éléments
        for idx, (i, j) in enumerate(self.connectivity):
            elem_type = self.elem_types[idx]
            color = self.color_map.get(elem_type, 'black')
            
            # Structure initiale (pointillés)
            self.ax.plot(
                [self.x0[i], self.x0[j]], [self.y0[i], self.y0[j]], [self.z0[i], self.z0[j]],
                linestyle='--', color='gray', linewidth=0.5, alpha=0.5
            )
            
            # Structure déformée
            self.ax.plot(
                [xd[i], xd[j]], [yd[i], yd[j]], [zd[i], zd[j]],
                color=color, linewidth=2
            
            )

        # --- MODIFICATION ICI ---
        # Affichage des numéros de nœuds ET déplacements
        if self.show_nodes:
            # Conversion de initial_point_index en liste pour une recherche rapide "in"
            target_indices = np.atleast_1d(self.initial_point_index).tolist()
            
            for n in range(len(self.x0)):
                # Texte de base : le numéro du nœud
                label_text = f'{n+1}'
                text_color = 'blue'
                
                # Si le nœud fait partie des points initiaux demandés
                if n in target_indices:
                    # On récupère les valeurs réelles (sans scale)
                    dx = self.Ux[n]
                    dy = self.Uy[n]
                    dz = self.Uz[n]
                    # On ajoute les déplacements au texte (notation scientifique pour la lisibilité)
                    label_text += f"\nUx:{dx:.2e}\nUy:{dy:.2e}\nUz:{dz:.2e}"
                    text_color = 'red' # On change la couleur pour les mettre en évidence

                self.ax.text(xd[n], yd[n], zd[n], label_text, color=text_color, fontsize=8)
        # ------------------------

        # Mise à l'échelle des axes (Equal Aspect Ratio)
        all_x = np.concatenate([self.x0, xd])
        all_y = np.concatenate([self.y0, yd])
        all_z = np.concatenate([self.z0, zd])
        
        max_range = np.array([all_x.max()-all_x.min(), 
                              all_y.max()-all_y.min(), 
                              all_z.max()-all_z.min()]).max() / 2.0

        mid_x = (all_x.max()+all_x.min()) * 0.5
        mid_y = (all_y.max()+all_y.min()) * 0.5
        mid_z = (all_z.max()+all_z.min()) * 0.5

        self.ax.set_xlim(mid_x - max_range, mid_x + max_range)
        self.ax.set_ylim(mid_y - max_range, mid_y + max_range)
        self.ax.set_zlim(mid_z - max_range, mid_z + max_range)

        # Labels et Titres
        self.ax.set_xlabel('X [m]')
        self.ax.set_ylabel('Y [m]')
        self.ax.set_zlabel('Z [m]')
        self.ax.set_title("Visualisation de la deformation")
      
        # Texte d'info (Scale et Lambda)
        if self.convergence:
            info_text = f"Convergence"
        else:
            info_text = f"Non-convergence"

        self.ax.text2D(0.05, 0.95, info_text, transform=self.ax.transAxes, color='black')

        self.fig.canvas.draw_idle()


def post_process(Coord, U, Connect, Elem_Types, initial_point_index, _lambda):
    """
    Fonction wrapper pour conserver la signature originale.
    """
    # Instanciation de la classe de visualisation
    viz = StructureVisualizer(Coord, U, Connect, Elem_Types, initial_point_index, _lambda)
    plt.show()