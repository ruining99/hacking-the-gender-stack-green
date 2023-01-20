import { Outlet, Route, createRoutesFromElements } from 'react-router';

import ErrorPage from './error-page/ErrorPage';
import Grammarly from './grammarly/grammarly.jsx';
import Home from './home/Home';
import { withMainLayout } from './layout/main/Main';


export default createRoutesFromElements(
  <Route path="/" element={withMainLayout(<Outlet />)} errorElement={withMainLayout(<ErrorPage />)}>
    <Route index element={<Home />} />
<<<<<<< HEAD
    <Route index path='grammarly' element={<Grammarly />} />
=======
    <Route path="grammarly" element={<div />} />
>>>>>>> cabd3dea6eff45b6cb660181dea7b2e2ece5ddf0
  </Route>
);
