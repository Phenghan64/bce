"""Inference-only adapters for the published AnOxPePred and EmulsiPred models.

No fitting or score calibration. See upstream archives and SHA256 manifest.
"""
from pathlib import Path
import importlib.util
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ANOX = ROOT / '.tools/enzyme-case/AnOxPePred-master'
EMUL = ROOT / '.tools/enzyme-case/EmulsiPred-master/EmulsiPred'


class PublishedModels:
    def __init__(self):
        import tensorflow as tf
        tf.config.threading.set_intra_op_parallelism_threads(4)
        tf.config.threading.set_inter_op_parallelism_threads(2)
        class AnOxPePred(tf.keras.Model):
            def __init__(self):
                super().__init__()
                self.conv = tf.keras.layers.Conv1D(128, 3, strides=1, activation='elu', padding='same')
                self.maxpool = tf.keras.layers.AveragePooling1D(pool_size=3, strides=3)
                self.dropout1 = tf.keras.layers.Dropout(.1)
                self.flatten = tf.keras.layers.Flatten()
                self.d1 = tf.keras.layers.Dense(256, activation='elu')
                self.dropout2 = tf.keras.layers.Dropout(.15)
                self.d2 = tf.keras.layers.Dense(2, activation='sigmoid')
            def call(self, x):
                return self.d2(self.dropout2(self.d1(self.flatten(self.dropout1(self.maxpool(self.conv(x)), training=False))), training=False))
        self.model = AnOxPePred()
        self.model(np.zeros((1, 30, 20), dtype=np.float32))
        # Explicitly assign all six tensors: old TF1 object configs are not TF2 variables.
        checkpoint = tf.train.load_checkpoint(str(ANOX / 'Result/AnOxPePred_v1'))
        for name in ('conv', 'd1', 'd2'):
            weights = [checkpoint.get_tensor(f'{name}/{v}/.ATTRIBUTES/VARIABLE_VALUE') for v in ('kernel', 'bias')]
            layer = getattr(self.model, name)
            assert [tuple(x.shape) for x in weights] == [tuple(x.shape) for x in layer.weights]
            layer.set_weights(weights)
        self.embedding = {}
        for line in (ANOX / 'Data/One-hot_encoding.txt').read_text().splitlines():
            if line and not line.startswith('#'):
                x = line.split(); self.embedding[x[0]] = list(map(float, x[1:]))
        spec = importlib.util.spec_from_file_location('original_emulsi_utils', EMUL / 'utils.py')
        self.emul = importlib.util.module_from_spec(spec); spec.loader.exec_module(self.emul)
        self.hydro = self.emul.hydrophobicity_scale(1)
        self.norm = {}
        for k in ('a', 'b', 'g'):
            frame = pd.read_csv(EMUL / f'NormalizationValues/{k}_norm.csv', index_col=0)
            self.norm[k] = self.emul.gamma_norm_df_to_norm_dict(frame) if k == 'g' else {int(n): (m, s) for n,m,s in frame.values}

    def anox(self, sequences, batch_size=256):
        output = []
        for begin in range(0, len(sequences), batch_size):
            batch = sequences[begin:begin+batch_size]
            assert all(3 <= len(s) <= 30 and set(s) <= set('ACDEFGHIKLMNPQRSTVWY') for s in batch)
            padded = ['X'*((30-len(s)+1)//2)+s+'X'*((30-len(s))//2) for s in batch]
            x = np.array([[self.embedding[c] for c in s] for s in padded], dtype=np.float32)
            output.extend(self.model(x, training=False).numpy().tolist())
        return np.array(output)

    def emulsify(self, s):
        assert 7 <= len(s) <= 30
        return [float(self.emul.calculate_emul(s, self.hydro, self.emul.alpha_emul, self.norm['a'])),
                float(self.emul.calculate_emul(s, self.hydro, self.emul.beta_emul, self.norm['b'])),
                float(self.emul.calculate_gamma_emul(s, self.hydro, self.norm['g'])[0])]


if __name__ == '__main__':
    import json
    m = PublishedModels()
    anchors = {'IIAPPERKYS': (.40,.23), 'GADPEDVIVA': (.32,.20), 'DIDIRKDLYAN': (.27,.17), 'VAPEEHPTL': (.51,.25), 'LLPVLYPPVVEE': (.64,.21)}
    out = []
    for s, pred in zip(anchors, m.anox(list(anchors))):
        em = m.emulsify(s)
        out.append({'sequence':s,'anox':pred.tolist(),'paper_anox':anchors[s],'emulsi':em})
        assert np.max(np.abs(pred - np.array(anchors[s]))) < .011, (s, pred)
    assert abs(out[0]['emulsi'][2]-3.43)<.02
    assert abs(out[2]['emulsi'][1]-3.65)<.02
    (ROOT/'results/enzyme-case/model-compatibility.json').write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
