import {inject} from '@angular/core';
import { ActivatedRouteSnapshot, RouterStateSnapshot } from "@angular/router";
import {Observable} from "rxjs";
import {MaveDBDetailsResponse} from "../model/mavedb";
import {ApiService} from "./api.service";


export function maveDBRecordResolver(
  route: ActivatedRouteSnapshot,
  state: RouterStateSnapshot): Observable<MaveDBDetailsResponse> {

  const apiService = inject(ApiService);

  return apiService.getMaveDBRecord(route.params['urn']);

}
