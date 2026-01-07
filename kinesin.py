import numpy as np

class Kinesin:
    
    class KHead:
        def __init__(self, size, state, start, D, dt, kb, T, s, k_atp_bind, k_atp_hyd, k_adp_release):
            self.pos = np.zeros(size)
            self.rand = np.sqrt(2*D*dt) * np.random.normal(loc=0, scale=1, size=size)
            self.state = state
            self.pos[0] = start
            self.V0 = kb * T
            self.s = s
            self.k = (dt * D) / (kb * T) # dt/g
            self.k_atp_bind = k_atp_bind
            self.k_adp_release = k_adp_release
            self.k_atp_hyd = k_atp_hyd
            self.state_history = np.zeros(size)
            self.state_history[0] = self.state_index()
            self.schedule_next_state(0)
            
        def state_index(self):
            match(self.state):
                case 'T':
                    return 2
                case 'D':
                    return 3
                case 'E':
                    return 1

        def binding_force(self, counter):
            n = int(self.pos[counter-1]/8)
            binding_sites = np.array([i*8 for i in range(n-4, n+5)])
            X = self.pos[counter-1] - binding_sites
            exp_sum = 2 * np.sum(np.exp(-X**2) * X)

            match(self.state):
                case 'T':
                    mVb0 = 16
                case 'D':
                    mVb0 = 7
                case 'E':
                    mVb0 = 10

            return self.V0 * (-mVb0 * exp_sum + self.s)

        def attractive_force(self, other, counter):
            max_separation = 9
            separation = self.pos[counter-1] - other.pos[counter-1]
            abs_sep = np.abs(separation)
            if abs_sep == 0:
                return 0
            t1 = abs_sep - max_separation
            return -self.V0 * np.exp(t1) * separation / abs_sep

        def repulsive_force(self, other, counter):
            separation = self.pos[counter-1] - other.pos[counter-1]
            return 7 * 2 * self.V0 * np.exp(-separation**2) * separation

        def update_pos(self, other, counter):
            Fb = self.binding_force(counter)
            Fatt = self.attractive_force(other, counter)
            Frep = self.repulsive_force(other, counter)
            Ftotal = Fb + Fatt + Frep
            step = Ftotal * self.k + self.rand[counter]
            self.pos[counter] = self.pos[counter-1] + step

        def schedule_next_state(self, current_time):
            match(self.state):
                case 'T':
                    rate = self.k_atp_hyd
                    self.next_state = 'D'
                case 'D':
                    rate = self.k_adp_release
                    self.next_state = 'E'
                case 'E':
                    rate = self.k_atp_bind
                    self.next_state = 'T'

            self.time_next_event = current_time - np.log(np.random.rand()) / (rate + 1e-50)

        def update_state(self, other, current_time, counter):
            if current_time >= self.time_next_event:
                self.state = self.next_state
                self.schedule_next_state(current_time)
            self.state_history[counter] = self.state_index()

    def __init__(self, size, startl, startr, statel, stater, D, dt, kb, T, s, atp_conc, k_atp_bind, k_atp_hyd, k_adp_release):
        size += 1
        self.headl = self.KHead(size, statel, startl, D, dt, kb, T, s, atp_conc*k_atp_bind, k_atp_hyd, k_adp_release)
        self.headr = self.KHead(size, stater, startr, D, dt, kb, T, s, atp_conc*k_atp_bind, k_atp_hyd, k_adp_release)
        self.time = np.zeros(size)
        self.counter = 1
        self.size = size
        self.dt = dt

    def update_pos(self):
        self.headl.update_pos(self.headr, self.counter)
        self.headr.update_pos(self.headl, self.counter)
        self.time[self.counter] = self.time[self.counter-1] + self.dt

    def update_state(self):
        self.headl.update_state(self.headr, self.time[self.counter], self.counter)
        self.headr.update_state(self.headl, self.time[self.counter], self.counter)

    def update(self):
        self.update_pos()
        self.update_state()
        self.counter += 1

    def data(self):
        cutdown = int(self.size / 100000)
        self.time /= 1000 # to ms
        return [self.time[::cutdown], self.headl.pos[::cutdown], self.headr.pos[::cutdown], self.headl.state_history[::cutdown], self.headr.state_history[::cutdown]]
