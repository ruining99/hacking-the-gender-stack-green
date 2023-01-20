import checker from 'vite-plugin-checker';
import react from '@vitejs/plugin-react';
import fs from 'fs';
import module from 'module';
import path from 'path';

const WRONG_CODE = `import { bpfrpt_proptype_WindowScroller } from "../WindowScroller.js";`;

export function reactVirtualized() {
  return {
    name: 'flat:react-virtualized',
    // Note: we cannot use the `transform` hook here
    //       because libraries are pre-bundled in vite directly,
    //       plugins aren't able to hack that step currently.
    //       so instead we manually edit the file in node_modules.
    //       all we need is to find the timing before pre-bundling.
    configResolved() {
      const require = module.createRequire(import.meta.url);
      const file = require
        .resolve('react-virtualized')
        .replace(
          path.join('dist', 'commonjs', 'index.js'),
          path.join('dist', 'es', 'WindowScroller', 'utils', 'onScroll.js')
        );
      const code = fs.readFileSync(file, 'utf-8');
      const modified = code.replace(WRONG_CODE, '');
      fs.writeFileSync(file, modified);
    },
  };
}

/**
 * @type {import('vite').UserConfigExport}
 */
export default {
  plugins: [
    react(),
    checker({
      eslint: { lintCommand: 'eslint src/**/*.{ts,tsx}' },
      stylelint: { lintCommand: 'stylelint src/**/*.scss' },
      typescript: true,
    }),
    reactVirtualized(),
  ],
  server: {
    port: 3000,
    hmr: {
      /**
       * Port forwarding in codespaces doesn't open the local dev server port on the proxy so HMR
       * websocket connections should just be made via the HTTPS port
       */
      clientPort: process.env.CODESPACES ? 443 : 3000,
    },

    proxy: {
      '/api': 'http://127.0.0.1:8000',
    },
  },
};
